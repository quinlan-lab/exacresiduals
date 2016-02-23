# SEE: https://github.com/quinlan-lab/lab-wiki/blob/master/projects/residuals.md
import sys
import os
import gzip
from collections import defaultdict
import subprocess
import operator
import itertools as it
from bisect import bisect_left

from cyvcf2 import VCF
import numpy as np
from bw import BigWig

xopen = lambda f: (gzip.open if f.endswith(".gz") else open)(f)


def path(p):
    return os.path.expanduser(os.path.expandvars(p))


def read_gerp(chrom, gerp_path=path("~u6000771/Projects/gemini_install/data/gemini_data/hg19.gerp.bw")):
    gerp = BigWig(gerp_path)
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    d = dict(gerp.chroms)
    l = d[chrom]

    return np.frombuffer(gerp.values(chrom, 0, l), dtype='f')


def read_coverage(chrom, cov=10, length=249250621, path="~u6000771/Data/ExAC-coverage/"):
    """
    read ExAC coverage from a single chrom into a numpy array. If no length is
    given, just use the one length from chrom 1.
    path is expected to contain Panel.chr*
    cov is the column to pull
    """

    cols = "chrom	pos	mean	median	1	5	10	15	20	25	30	50 100".split()
    coli = cols.index(str(cov)) + 1

    # just extract the position (2) and the requested column
    p = subprocess.Popen("tabix {path}/Panel.chr{chrom}.coverage.txt.gz {chrom} | cut -f 2,{coli} ".format(**locals()),
            stdout=subprocess.PIPE, stderr=sys.stderr,
            shell=True,
            executable=os.environ.get("SHELL"))

    cov = np.zeros(length, dtype=np.float32)
    j = 0
    for line in p.stdout:
        pos, val = line.split()
        cov[int(pos)-1] = float(val)
        j += 1
    assert j > 0, ("no values found for", chrom, path)
    p.wait()
    if p.returncode != 0:
        raise Exception("bad: %d", p.returncode)
    return cov


def read_exons(gtf):
    transcripts = defaultdict(list)
    ends = defaultdict(list)
    for toks in (x.rstrip('\r\n').split("\t") for x in xopen(gtf) if x[0] != "#"):
        if toks[2] not in("UTR", "exon"): continue
        start, end = map(int, toks[3:5])
        transcript = toks[8].split('transcript_id "')[1].split('"', 1)[0]
        transcripts[transcript].append(start)
        ends[transcript].append(end)

    # sort by start so we can do binary search.
    transcripts = dict((k, sorted(v)) for k, v in transcripts.iteritems())
    return transcripts, ends


def get_cdna_start_end(cdna_start):
    cdna_start = cdna_start.rstrip("-?")
    if cdna_start[0] == "?": # deletion
        _, cdna_end = cdna_start.split("-")
        cdna_end = int(cdna_end)
        cdna_start = cdna_end - len(v.REF)
    elif "-" in cdna_start:
        try:
            cdna_start, cdna_end = map(int, cdna_start.split("-"))
        except:
            print v.REF, v.ALT, cdna_start, csq
            raise
    else:
        cdna_start = int(cdna_start)
        cdna_end = cdna_start + len(v.REF)
    return cdna_start, cdna_end

exac = VCF('/usr/local/src/gemini_install/data/gemini_data/ExAC.r0.3.sites.vep.tidy.vcf.gz')

# CSQ keys
kcsq = exac["CSQ"]["Description"].split(":")[1].strip(' "').split("|")


def isfunctional(csq):
    return any(c in ('stop_gained', 'stop_lost', 'start_lost', 'initiator_codon_variant', 'rare_amino_acid_variant',
                     'missense_variant', '5_prime_UTR_premature_start_codon_gain_variant', 'protein_altering_variant',
                     'splice_acceptor_variant', 'splice_donor_variant', 'frameshift_variant')
               for c in csq['Consequence'].split('&'))

# read ensembl gtf into dict keyed by transcript with list of exons so
# we know how far back to go.
transcript_exon_starts, transcript_exon_ends = read_exons("Homo_sapiens.GRCh37.75.gtf.gz")


header = "chrom\tstart\tend\taf\tfunctional\tgene\ttranscript\texon\timpact\tcdna_start\tcdna_end\tcoverage\tgerp\tposns"
print "#" + header
keys = header.split("\t")
for chrom, viter in it.groupby(exac, operator.attrgetter("CHROM")):
    rows = []

    #gerp_array = read_gerp(chrom) # handle nans?
    #coverage_array = read_coverage(chrom, depth=10)

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


        cdna_start, cdna_end = get_cdna_start_end(csq['cDNA_position'])

        rows.append(dict(chrom=v.CHROM, start=v.start, end=v.end, af=af,
            functional=int(isfunctional(csq)),
            gene=csq['SYMBOL'], transcript=csq['Feature'], exon=csq['EXON'],
            impact=csq['Consequence'],
            cdna_start=cdna_start,   cdna_end=cdna_end))

        # TODO: remove this. just testing on a subset.
        if len(rows) > 100000: break

    # now we need to sort and then group by transcript so we know the gaps.
    rows.sort(key=operator.itemgetter('transcript', 'start'))

    out = []
    for transcript, trows in it.groupby(rows, operator.itemgetter("transcript")):
        last = 0

        exon_starts = transcript_exon_starts[transcript]

        for i, row in enumerate(trows, start=1):


            istart = bisect_left(exon_starts, last)
            iend = bisect_left(exon_starts, row['cdna_start'])

            # easy case is when variants are in in same exon; just grab all the
            # scores with a single query
            diff = row['cdna_start'] - last
            if istart == iend:
                qstart, qend = row['start'] - diff, row['start']
                row['gerp'] = "XX" # ",".join("%.2f" % g for g in gerp_array[qstart:qend])
                row['coverage'] = "XX" # ",".join("%.2f" % g for g in coverage_array[qstart:qend])
                row['posns'] = ",".join(map(str, range(qstart, qend)))
            else:
                # loop over exons until we have queried diff bases.
                L_gerp, L_coverage, L_posns = [], [], []
                exon_ends = transcript_exon_ends[transcript]
                for xstart, xend in zip(exon_starts[istart-1:], exon_ends[istart-1:]):
                    assert xstart <= xend
                    # had to go to start of exon so we take the max but this is
                    # only required for the 1st time through the loop.
                    xstart = max(xstart, last)
                    xend = min(xend, xstart + diff) # dont read more than we need
                    #L_gerp.extend("%.2f" % g for g in gerp_array[xstart:xend])
                    #L_coverage.extend("%.2f" % g for g in coverage_array[xstart:xend])
                    L_posns.extend(str(s) for s in range(xstart, xend))

                    if len(L_posns) >= diff: break

            # TODO:
            # when i == len(trows) add an extra column to get to end of
            # transcript? or extra row?

            last = row['cdna_end'] # or start?

            out.append(row)

    # still print in sorted order
    out.sort(key=operator.itemgetter('start'))
    for d in out:
        print "\t".join(map(str, (d[k] for k in keys)))
