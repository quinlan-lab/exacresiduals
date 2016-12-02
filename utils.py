import sys
import os
from collections import defaultdict
from bisect import bisect_left
import subprocess
import toolshed as ts

from interlap import InterLap, Interval as IntervalSet, reduce as ireduce
import numpy as np

def split_ranges(position, ranges, splitters): # if range is in splitters, it is removed from potential constraint regions
    """
    >>> split_ranges(1033, [(1018, 1034)], [(1022, 1034)])
    [[(1018, 1022)]]

    >>> split_ranges(1033, [(1018, 1034)], None)
    [[(1018, 1034)]]

    >>> split_ranges(1033, [(1018, 1034)], [(1022, 1024), (1028, 1034)])
    [[(1018, 1022)], [(1024, 1028)]]

    >>> split_ranges(57, [(18, 24), (28, 35), (55, 60)], [(28, 35), (55, 57)])
    [[(18, 24)], [(57, 60)]]

    >>> split_ranges(5, [(12, 18), (22, 28), (32, 39), (42, 48)],
    ...                 [(12, 18), (22, 26),           (44, 48)])
    [[(26, 28)], [(32, 39)], [(42, 44)]]

    >>> split_ranges(5, [(11, 18), (22, 28), (32, 39), (42, 48)],
    ...                 [(12, 18), (22, 26),           (44, 48)])
    [[(11, 12)], [(26, 28)], [(32, 39)], [(42, 44)]]

    """
    if splitters is None:
        return [ranges]

    return [x._vals for x in IntervalSet(ranges).split(splitters)]

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
    assert last+1 >= exon_starts[0]
    assert vstart <= exon_ends[-1]
    assert vstart+1 >= last, (vstart, last, exon_starts)
    assert all(s < e for s, e in zip(exon_starts, exon_ends))

    istart = bisect_left(exon_starts, last) - 1
    if istart == -1: istart = 0
    [(0, 6), (10, 11)]

    assert exon_starts[istart] <= last+1, (exon_starts[istart], last, istart)
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

import doctest
res = doctest.testmod()
if res.failed != 0:
    sys.stderr.write("FAILING TESTS")
    sys.exit(1)


def path(p):
    return os.path.expanduser(os.path.expandvars(p))

def floatfmt(v, prec="%.2f"):
    return (prec % v).rstrip('0').rstrip('.')

def read_coverage(chrom, cov=10, length=249250621, path="data/"):
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


def read_exons(gtf, chrom, coverage_array, *args):
    genes = defaultdict(IntervalSet)
    splitters = defaultdict(IntervalSet)


    interlaps = []
    split_iv = InterLap()
    # preempt any bugs by checking that we are getting a particular chrom
    assert gtf[0] == "|", ("expecting a tabix query so we can handle chroms correctly")
    f1 = open("selfchainremoved.txt","a")
    f2 = open("segdupsremoved.txt","a")
    f3 = open("coveragecut.txt","a")
    for a in args:
        assert a[0] == "|", ("expecting a tabix query so we can handle chroms correctly", a)
    
        # any file that gets sent in will be used to split regions (just like
        # low-coverage). For example, we split on self-chains as well.
        for toks in (x.strip().split("\t") for x in ts.nopen(a)): # adds self chains and segdups to splitters list, so that exons can be split, and they are removed from CCRs
            s, e = int(toks[1]), int(toks[2])
            split_iv.add((s, e))
            if len(toks) > 3:
                f1.write("\t".join(toks)+"\n")
            else:
                f2.write("\t".join(toks)+"\n")
                

    for toks in (x.rstrip('\r\n').split("\t") for x in ts.nopen(gtf) if x[0] != "#"):
        if toks[2] not in("CDS", "stop_codon") or toks[1] not in("protein_coding"): continue
        #if toks[0] != "1": break
        start, end = map(int, toks[3:5])
        gene = toks[8].split('gene_name "')[1].split('"', 1)[0]
        assert start <= end, toks
        key = toks[0], gene

        cutoff = 0.3

        # find sections of exon under certain coverage.
        if coverage_array[start-1:end].min() < cutoff: # doesn't bother to run these operations if there is not one bp below the cutoff
            #splitters[key].add([(start - 1, end)]) #this takes out the whole exon for one section of poor coverage
            a = coverage_array[start - 1: end]
            #print str(start-1),end,a
            is_under, locs = False, [] # generates "locs" for each exon"
            if a[0] < cutoff:
                locs.append([start - 1])
                is_under = True # so you can initialize is_under
            for pos, v in enumerate(a[1:], start=start): #enumerates positions in the coverage array starting at the beginning of the exon
                if v < cutoff:
                    if not is_under:
                        is_under = True
                        locs.append([pos - 1]) #start, coverage is in bed format, so pos-1 is necessary, since splitters are open left and right side
                else:
                    if is_under:
                        is_under = False
                        locs[-1].append(pos) #end
            if is_under:
                locs[-1].append(end) # in this case would end splitter at the end of the exon
            splitters[key].add(map(tuple, locs))
            for i in locs:
                f3.write(chrom+"\t"+"\t".join(map(str,i))+"\n")

        for s, e in split_iv.find((start - 1, end)):
            splitters[key].add([(s, e)])

        genes[key].add([(start-1, end)])
    # sort by start so we can do binary search.
    genes = dict((k, sorted(v._vals)) for k, v in genes.iteritems())
    #ends = dict((k, sorted(v)) for k, v in ends.iteritems())
    splits, starts, ends = {}, {}, {}
    splitters = dict(splitters)
    for chrom_gene, sends in genes.iteritems():
        starts[chrom_gene] = [s[0] for s in sends]
        ends[chrom_gene] = [s[1] for s in sends]
        if chrom_gene in splitters:
            splits[chrom_gene] = splitters[chrom_gene]._vals

    return starts, ends, splits


def get_cdna_start_end(cdna_start, v):
    cdna_start = cdna_start.rstrip("-?")
    if cdna_start[0] == "?": # deletion
        _, cdna_end = cdna_start.split("-")
        cdna_end = int(cdna_end)
        cdna_start = cdna_end - len(v.REF)
    elif "-" in cdna_start:
        try:
            cdna_start, cdna_end = map(int, cdna_start.split("-"))
        except:
            print(v.REF, v.ALT, cdna_start, csq)
            raise
    else:
        cdna_start = int(cdna_start)
        cdna_end = cdna_start + len(v.REF)
    return cdna_start, cdna_end

def isfunctional(csq):
    return any(c in ('stop_gained', 'stop_lost', 'start_lost', 'initiator_codon_variant', 'rare_amino_acid_variant',
                     'missense_variant', 'protein_altering_variant', 'frameshift_variant')
               for c in csq['Consequence'].split('&'))


def cg_content(seq):
    if len(seq) == 0: return 0.0
    return 2.0 * seq.count('CG') / len(seq)
