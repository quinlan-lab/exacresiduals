# SEE: https://github.com/quinlan-lab/lab-wiki/blob/master/projects/residuals.md
import sys
import os
import gzip
from collections import defaultdict
import subprocess
import operator
import itertools as it
from bisect import bisect_left
from itertools import chain

from cyvcf2 import VCF
import numpy as np
from bw import BigWig
from pyfaidx import Fasta
import pyinter

from bisect import bisect_left

def get_ranges(last, vstart, exon_starts, exon_ends):
    """
    >>> get_ranges(61018, 62029, (
    ... 60174, 60370, 60665, 60925, 62029, 62216, 62453,
    ... 62675, 63052, 63398, 63652, 63868, 64512, 64764,
    ... 65018, 65671), (60281,
    ... 60565, 60808, 61033, 62134, 62379, 62587, 62824,
    ... 63209, 63559, 63779, 64102, 64691, 64946, 65084,
    ... 65985))
    [(61018, 61034), (62029, 62030)]

    >>> get_ranges(56, 95, range(0, 1000, 10), range(5, 1000, 10))
    [(60, 66), (70, 76), (80, 86), (90, 96)]

    >>> get_ranges(1, 10, range(0, 100, 10), range(5, 100, 10))
    [(1, 6), (10, 11)]

    >>> get_ranges(0, 10, range(0, 100, 10), range(5, 100, 10))
    [(0, 6), (10, 11)]

    >>> get_ranges(50, 60, (10,),(100,))
    [(50, 61)]

    >>> get_ranges(50, 60, (10,),(100,))
    [(50, 61)]

    >>> get_ranges(1562576, 1562675,
    ... (1560665, 1560925, 1562029, 1562216, 1562453, 1562675, 1563052, 1563398, 1563652, 1563868, 1564512, 1564764, 1565018, 1565671),
    ... (1560808, 1561033, 1562134, 1562379, 1562587, 1562824, 1563209, 1563559, 1563779, 1564102, 1564691, 1564946, 1565084, 1565985))
    [(1562576, 1562588), (1562675, 1562676)]
    """

    assert last >= exon_starts[0]
    assert vstart <= exon_ends[-1]
    assert vstart >= last
    assert all(s < e for s, e in zip(exon_starts, exon_ends))

    istart = bisect_left(exon_starts, last) - 1
    if istart == -1: istart = 0
    [(0, 6), (10, 11)]

    assert exon_starts[istart] <= last, (exon_starts[istart], last, istart)
    if exon_ends[istart] <= last:
        istart += 1
        last = exon_starts[istart]

    start = last
    ranges = []
    while start <= vstart and istart < len(exon_starts):
        ranges.append((start, exon_ends[istart] + 1))
        istart += 1
        if ranges[-1][1] > vstart:
            ranges[-1] = (ranges[-1][0], vstart + 1)
            break
        start = exon_starts[istart]

    return ranges

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


def floatfmt(v, prec="%.2f"):
    return (prec % v).rstrip('0').rstrip('.')

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
        #if j > 100000: break
    assert j > 0, ("no values found for", chrom, path)
    p.wait()
    if p.returncode != 0:
        raise Exception("bad: %d", p.returncode)
    return cov


def read_exons(gtf):
    transcripts = defaultdict(pyinter.IntervalSet)

    for toks in (x.rstrip('\r\n').split("\t") for x in xopen(gtf) if x[0] != "#"):
        if toks[2] not in("UTR", "exon"): continue
        #if toks[0] != "1": break
        start, end = map(int, toks[3:5])
        assert start <= end, toks
        transcript = toks[8].split('transcript_id "')[1].split('"', 1)[0]
        transcripts[transcript].add(pyinter.closedopen(start-1, end))

    # sort by start so we can do binary search.
    # TODO: need to remove overlapping exons so we don't double-count
    transcripts = dict((k, sorted(v)) for k, v in transcripts.iteritems())
    #ends = dict((k, sorted(v)) for k, v in ends.iteritems())
    starts, ends = {}, {}
    for tr, ivset in transcripts.iteritems():
        sends = sorted(list(ivset))
        starts[tr] = [x.lower_value for x in sends]
        ends[tr] = [x.upper_value for x in sends]
    return starts, ends


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

exac = VCF('/uufs/chpc.utah.edu/common/home/u6000771/Projects/gemini_install/data/gemini_data/ExAC.r0.3.sites.vep.tidy.vcf.gz')
#exac = VCF('rbp7.vcf.gz')

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

fasta = Fasta('/uufs/chpc.utah.edu/common/home/u6000771/Data/data/hs37d5.fa', read_ahead=10000, as_raw=True)
def cg_content(seq):
    if len(seq) == 0: return 0.0
    return 2.0 * seq.count('CG') / len(seq)

header = "chrom\tstart\tend\taf\tfunctional\tgene\ttranscript\texon\timpact\tvstart\tvend\tn_bases\tcg_content\tcdna_start\tcdna_end\tranges\tcoverage\tgerp\tposns"
print "#" + header
keys = header.split("\t")
for chrom, viter in it.groupby(exac, operator.attrgetter("CHROM")):
    rows = []
    print >>sys.stderr, "reading chrom",

    fa = fasta[chrom]

    gerp_array = read_gerp(chrom)
    print >>sys.stderr, "gerp",
    coverage_array = read_coverage(chrom, length=len(gerp_array), cov=10)
    print >>sys.stderr, chrom

    for v in viter:


        if not (v.FILTER is None or v.FILTER == "PASS"):
            continue
        info = v.INFO
        try:
            csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]
        except KeyError:
            continue
        af = info['AC_Adj'] / float(info['AN_Adj'] or 1)
        for csq in (c for c in csqs if c['CANONICAL'] == 'YES' and c['Allele'] == v.ALT[0]):
            # skipping intronic
            if csq['Feature'] == '' or csq['EXON'] == '' or csq['cDNA_position'] == '': continue
            if not isfunctional(csq): continue

            cdna_start, cdna_end = get_cdna_start_end(csq['cDNA_position'])

            rows.append(dict(chrom=v.CHROM, vstart=v.start, vend=v.end, af=af,
                functional=int(isfunctional(csq)),
                gene=csq['SYMBOL'], transcript=csq['Feature'], exon=csq['EXON'],
                impact=csq['Consequence'],
                cdna_start=cdna_start,   cdna_end=cdna_end))


    # now we need to sort and then group by transcript so we know the gaps.
    rows.sort(key=operator.itemgetter('transcript', 'vstart', 'vend'))

    out = []
    for transcript, trows in it.groupby(rows, operator.itemgetter("transcript")):
        exon_starts = transcript_exon_starts[transcript]
        exon_ends = transcript_exon_ends[transcript]
        last = exon_starts[0]
        for i, row in enumerate(trows, start=1):
            # istart and iend determin if we need to span exons.

            assert row['vstart'] <= exon_ends[-1], (row, exon_ends)
            ranges = get_ranges(last, row['vstart'], exon_starts, exon_ends)

            row['gerp'] = ",".join(",".join(floatfmt(g) for g in gerp_array[s:e]) for s, e in ranges)
            row['coverage'] = ",".join(",".join(floatfmt(g) for g in coverage_array[s:e]) for s, e in ranges)
            row['posns'] = list(chain.from_iterable([range(s, e) for s, e in ranges]))
            row['ranges'] = ["%d-%d" % (s, e) for s, e in ranges]
            seqs = [fa[s-1:e] for s, e in ranges]
            # this can happend for UTR variants since we can't really get
            # anything upstream of them.
            if row['posns'] == []:  # UTR:
                p = row['vstart']
                row['gerp'] = ",".join(floatfmt(g) for g in gerp_array[p:p+1])
                row['coverage'] = ",".join(floatfmt(g) for g in coverage_array[p:p+1])
                row['posns'] = [p]

            # post-hoc sanity check
            exon_bases = set(chain.from_iterable(range(s, e) for s, e in zip(exon_starts, exon_ends)))
            ranges = set(chain.from_iterable(range(int(x[0]), int(x[1])) for x in (z.split("-") for z in row['ranges'])))
            m = len(ranges - exon_bases)
            if m > len(row['ranges']):
                print >>sys.stderr, last, row['vstart'], row['ranges'], len(ranges - exon_bases), zip(exon_starts, exon_ends)

            # start or end? if we use end then can have - diff.
            row['ranges'] = ",".join(row['ranges'])
            last = row['vstart']
            row['n_bases'] = len(row['posns'])
            row['start'] = str(min(row['posns']))
            row['end'] = str(max(row['posns']))
            row['posns'] = ",".join(map(str, row['posns']))
            row['cg_content'] = floatfmt(np.mean([cg_content(s) for s in seqs]))
            if row['cg_content'] == 'nan':
                row['cg_content'] = '0'

            out.append(row)

    # still print in sorted order
    out.sort(key=operator.itemgetter('start'))
    for d in out:
        print "\t".join(map(str, (d[k] for k in keys)))
