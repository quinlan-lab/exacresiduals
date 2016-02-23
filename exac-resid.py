# SEE: https://github.com/quinlan-lab/lab-wiki/blob/master/projects/residuals.md
import sys
from collections import defaultdict
import operator
import itertools as it

from cyvcf2 import VCF
exac = VCF('/usr/local/src/gemini_install/data/gemini_data/ExAC.r0.3.sites.vep.tidy.vcf.gz')

# CSQ keys
kcsq = exac["CSQ"]["Description"].split(":")[1].strip(' "').split("|")

exons = defaultdict(list)

def read_gerp(chrom):
    # TODO bigwig read

def read_coverage(chrom, depth=10):
    # TODO:


def isfunctional(csq):
    return any(c in ('stop_gained', 'stop_lost', 'start_lost', 'initiator_codon_variant', 'rare_amino_acid_variant',
                     'missense_variant', '5_prime_UTR_premature_start_codon_gain_variant', 'protein_altering_variant',
                     'splice_acceptor_variant', 'splice_donor_variant', 'frameshift_variant')
               for c in csq['Consequence'].split('&'))

# TODO: add a column that has all coverage values (% @10X for max(current.transcript.start, previous.end + 1):current.start
# TODO:                     ... same for GERP scores
# A. groupby chromosome, then sort by (transcript, start).
# B. read GERP chrom and coverage chrom into memory. then gerp = GERP[last_end:current_start]
header = "chrom\tstart\tend\taf\tfunctional\tgene\ttranscript\texon\timpact\tcnda_start\tcdna_end\tcoverage\tgerp"
print "#" + header
keys = header.split("\t")
for chrom, viter in it.groupby(exac, operator.attrgetter("CHROM")):
    rows = []

    gerp_array = read_gerp(chrom) # handle nans?
    coverage_array = read_coverage(chrom, depth=10)

    for v in viter:
        if not (v.FILTER is None or v.FILTER == "PASS"):
            continue
        info = v.INFO
        try:
            csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]
        except KeyError:
            continue
        af = info['AC_Adj'] / float(info['AN_Adj'] or 1)
        try:
            csq = next(c for c in csqs if c['CANONICAL'] == 'YES' and c['Allele'] == v.ALT[0])
        except StopIteration:
            continue

        # skipping intronic
        if csq['Feature'] == '' or csq['EXON'] == '': continue

        cdna_start = csq['cDNA_position']
        if "-" in cdna_start:
            cdna_start, cdna_end = map(int, cdna_start.split("-"))
        else:
            cdna_start = int(cdna_start)
            cdna_end = cdna_start + len(v.REF)

        rows.append(dict(chrom=v.CHROM, start=v.start, end=v.end, af=af,
            functional=int(isfunctional(csq)),
            gene=csq['SYMBOL'], transcript=csq['Feature'], exon=csq['EXON'],
            impact=csq['Consequence'],
            cdna_start=cdna_start,   cdna_end=cdna_end))

    # now we need to sort and then group by transcript so we know the gaps.
    rows.sort(key=operator.itemgetter('transcript', 'start'))

    out = []
    for transcript, trows in it.groupby(rows, operator.itemgetter("transcript")):
        last = 0
        for i, row in enumerate(trows, start=1):
            diff = row['cdna_start'] - last
            qstart, qend = row['start'] - diff, row['start']

            row['gerp'] = ",".join("%.2f" % g for g in gerp_array[qstart:qend])
            row['coverage'] = ",".join("%.2f" % g for g in coverage_array[qstart:qend])

            # TODO:
            # when i == len(trows) add an extra column to get to end of
            # transcript? or extra row?

            last = row['cdna_end']  # or start?

            out.append(row)

    # still print in sorted order
    out.sort(key=operator.itemgetter('start'))
    for d in out:
        print "\t".join(map(str, (d[k] for k in keys)))
