from __future__ import print_function

# ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/ExAC.r0.3.sites.vep.vcf.gz
VCF_PATH = "data/ExAC.r0.3.sites.vep.vcf.gz"

# ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
GTF_PATH = "data/Homo_sapiens.GRCh37.75.gtf.gz"

# ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/coverage
COVERAGE_PATH = "data/"

# from UCSC. see data/get-chain.py, pipe output to sort -k1,1 -k2,2n | uniq | bgzip -c > data/self-chains.gt90.bed.gz
SELF_CHAINS = "data/self-chains.gt90.bed.gz"

# from UCSC.  data/segmental.bed.gz

SEGDUPS = "data/segmental.bed.gz"

FASTA_PATH = "/uufs/chpc.utah.edu/common/home/u6000771/Data/data/hs37d5.fa"


import sys
import itertools as it
import operator

import numpy as np
import utils as u
from cyvcf2 import VCF
from pyfaidx import Fasta

zip = it.izip


exac = VCF(VCF_PATH)
kcsq = exac["CSQ"]["Description"].split(":")[1].strip(' "').split("|")

#exac = exac("2:112538945-112551053")


fasta = Fasta(FASTA_PATH, read_ahead=10000, as_raw=True)

header = "chrom\tstart\tend\taf\tfunctional\tgene\ttranscript\texon\timpact\tvstart\tvend\tn_bases\tcg_content\tcdna_start\tcdna_end\tranges\tcoverage\tposns"
print("#" + header)
keys = header.split("\t")
for chrom, viter in it.groupby(exac, operator.attrgetter("CHROM")):
    rows = []
    print("reading chrom", file=sys.stderr)

    fa = fasta[chrom]
    coverage_array = u.read_coverage(chrom, length=len(fa), cov=10,
            path=COVERAGE_PATH)

    gene_exon_starts, gene_exon_ends, splitters = u.read_exons("|tabix {gtf} {chrom}"
                                                            .format(chrom=chrom,
                                                                gtf=GTF_PATH),
                                                            coverage_array,
                                            "|tabix {bed} {chrom}".format(chrom=chrom, bed=SELF_CHAINS),"|tabix {bed} {chrom}".format(chrom=chrom, bed=SEGDUPS))

    print(chrom, file=sys.stderr)
    for v in viter:
        if not (v.FILTER is None or v.FILTER == "PASS"):
            continue
        info = v.INFO
        is_multi = bool(info.get('OLD_MULTIALLELIC'))
        try:
            csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]
        except KeyError:
            continue
        # NOTE: using max here for alternates to be conservative
        ac = info['AC_Adj']
        if not isinstance(ac, (int, long)):
            ac = max(ac)
        af = ac / float(info['AN_Adj'] or 1)
        # NOTE: not requiring canonical or requiring the csq to match the
        # particular alt that we chose.
        for csq in (c for c in csqs if c['BIOTYPE'] == 'protein_coding'):
            # skipping intronic
            if csq['Feature'] == '' or csq['EXON'] == '' or csq['cDNA_position'] == '': continue
            if not u.isfunctional(csq): continue

            cdna_start, cdna_end = u.get_cdna_start_end(csq['cDNA_position'], v)

            rows.append(dict(chrom=v.CHROM, vstart=v.start, vend=v.end, af=af,
                functional=int(u.isfunctional(csq)),
                gene=csq['SYMBOL'], transcript=csq['Feature'], exon=csq['EXON'],
                impact=csq['Consequence'],
                cdna_start=cdna_start,   cdna_end=cdna_end))

    # now we need to sort and then group by gene so we know the gaps.
    rows.sort(key=operator.itemgetter('gene', 'vstart', 'vend'))

    out = []
    for chrom_gene, trows in it.groupby(rows, lambda row: (row['chrom'], row['gene'])):
        exon_starts = gene_exon_starts[chrom_gene]
        exon_ends = gene_exon_ends[chrom_gene]
        last = exon_starts[0]

        splitter = splitters.get(chrom_gene, None)

        for i, row in enumerate(trows, start=1):
            # istart and iend determin if we need to span exons.

            assert row['vstart'] <= exon_ends[-1], (row, exon_ends)
            mranges = u.get_ranges(last, row['vstart'], exon_starts, exon_ends)

            for ranges in u.split_ranges(row['vstart'], mranges, splitter):

                row['coverage'] = ",".join(",".join(u.floatfmt(g) for g in coverage_array[s:e]) for s, e in ranges)
                row['posns'] = list(it.chain.from_iterable([range(s, e) for s, e in ranges]))
                row['ranges'] = ["%d-%d" % (s, e) for s, e in ranges]
                seqs = [fa[s:e] for s, e in ranges]
                # this can happen for UTR variants since we can't really get
                # anything upstream of them.
                if row['posns'] == []:  # UTR:
                    p = row['vstart']
                    row['coverage'] = ",".join(u.floatfmt(g) for g in coverage_array[p:p+1])
                    row['posns'] = [p]

                # post-hoc sanity check
                exon_bases = set(it.chain.from_iterable(range(s, e) for s, e in zip(exon_starts, exon_ends)))
                ranges = set(it.chain.from_iterable(range(int(x[0]), int(x[1])) for x in (z.split("-") for z in row['ranges'])))
                m = len(ranges - exon_bases)
                if m > len(row['ranges']):
                    print(last, row['vstart'], row['ranges'], len(ranges -
                        exon_bases), zip(exon_starts, exon_ends),
                        file=sys.stderr)

                # start or end? if we use end then can have - diff.
                row['ranges'] = ",".join(row['ranges'])
                row['n_bases'] = len(row['posns'])
                row['start'] = str(min(row['posns']))
                row['end'] = str(max(row['posns']))
                row['posns'] = ",".join(map(str, row['posns']))
                row['cg_content'] = u.floatfmt(np.mean([u.cg_content(s) for s in seqs]))
                if row['cg_content'] == 'nan':
                    row['cg_content'] = '0'

                # we are re-using the dict for each loop so force a copy.
                out.append(dict(row))
            last = row['vstart']

    # still print in sorted order
    out.sort(key=operator.itemgetter('start'))
    last = (None, None)
    for d in out:
        key = d['start'], d['end']
        if key == last:
            continue
        last = key
        print("\t".join(map(str, (d[k] for k in keys))))
