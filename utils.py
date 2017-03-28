import sys
import os
from collections import defaultdict
from bisect import bisect_left
import subprocess
import toolshed as ts

from interlap import InterLap, Interval as IntervalSet, reduce as ireduce
import numpy as np


def overlaps(s1, e1, s2, e2):
    """
    >>> overlaps(2, 4, 3, 5)
    True
    >>> overlaps(2, 4, 4, 5)
    False
    >>> overlaps(2, 200, 3, 5)
    True
    >>> overlaps(3, 5, 2, 200)
    True
    >>> overlaps (3, 5, 5, 6)
    False
    >>> overlaps (5, 6, 3, 5)
    False
    """
    return not (e1 <= s2 or s1 >= e2)

def split_ranges(ranges, splitters, varflags): # if range is in splitters, it is removed from potential constraint regions; import my version of interlap
    """
    >>> split_ranges([(1018, 1034)], [(1022, 1034)], ['VARFALSE'])
    ([[(1018, 1022)]], ['VARFALSE'])

    >>> split_ranges([(1018, 1034)], [(1030, 1032)], ['VARFALSE'])
    ([[(1018, 1030)], [(1032, 1034)]], ['VARFALSE', 'VARFALSE'])

    >>> split_ranges([(1018, 1034)], None, ['VARFALSE'])
    ([[(1018, 1034)]], ['VARFALSE'])

    >>> split_ranges([(1018, 1034), (1045, 1069)], None, ['VARFALSE', 'VARTRUE'])
    ([[(1018, 1034)], [(1045, 1069)]], ['VARFALSE', 'VARTRUE'])

    >>> split_ranges([(1018, 1034), (1045, 1069)], [(1030, 1032), (1047, 1050)], ['VARFALSE', 'VARTRUE'])
    ([[(1018, 1030)], [(1032, 1034)], [(1045, 1047)], [(1050, 1069)]], ['VARFALSE', 'VARFALSE', 'VARTRUE', 'VARTRUE'])

    >>> split_ranges([(1018, 1034), (1045, 1069)], [(1030, 1050)], ['VARFALSE', 'VARTRUE'])
    ([[(1018, 1030)], [(1050, 1069)]], ['VARFALSE', 'VARTRUE'])

    >>> split_ranges([(1018, 1034)], [(1022, 1024), (1028, 1034)], ['VARFALSE'])
    ([[(1018, 1022)], [(1024, 1028)]], ['VARFALSE', 'VARFALSE'])

    >>> split_ranges([(18, 24), (28, 35), (55, 60)], [(28, 35), (55, 57)], ['VARFALSE', 'VARTRUE', 'VARTRUE'])
    ([[(18, 24)], [(57, 60)]], ['VARFALSE', 'VARTRUE'])

    >>> split_ranges([(12, 18), (22, 28), (32, 39), (42, 48)],
    ...                 [(12, 18), (22, 26),           (44, 48)], ['VARFALSE', 'VARFALSE', 'VARFALSE', 'VARTRUE'])
    ([[(26, 28)], [(32, 39)], [(42, 44)]], ['VARFALSE', 'VARFALSE', 'VARTRUE'])

    >>> split_ranges([(11, 18), (22, 28), (32, 39), (42, 48)],
    ...                 [(12, 18), (22, 26),           (44, 48)], ['VARFALSE', 'VARFALSE', 'VARFALSE', 'VARTRUE'])
    ([[(11, 12)], [(26, 28)], [(32, 39)], [(42, 44)]], ['VARFALSE', 'VARFALSE', 'VARFALSE', 'VARTRUE'])

    >>> split_ranges([(38826782, 38826890), (38827874, 38828144), (38828232, 38828286), (38834219, 38834405), (38834632, 38834759), (38834935, 38835008), (38837089, 38837266), (38842900, 38843135), (38845349, 38845487), (38847123, 38847183), (38847381, 38847501), (38848917, 38848968), (38849071, 38849201), (38850109, 38850221), (38851128, 38851287), (38851370, 38851483), (38852287, 38852501), (38852849, 38852912), (38853015, 38853214), (38853353, 38853472), (38855530, 38855611), (38855700, 38855757), (38857795, 38857952), (38858148, 38858212), (38858320, 38858416), (38858687, 38858777), (38860611, 38860612)],
    ...             [(38826782, 38826890), (38827874, 38828144), (38828232, 38828286), (38834219, 38834405), (38834632, 38834759), (38834935, 38835008), (38837089, 38837266), (38842900, 38843135)],
    ...             [ 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARTRUE'])
     ([[(38845349, 38845487)], [(38847123, 38847183)], [(38847381, 38847501)], [(38848917, 38848968)], [(38849071, 38849201)], [(38850109, 38850221)], [(38851128, 38851287)], [(38851370, 38851483)], [(38852287, 38852501)], [(38852849, 38852912)], [(38853015, 38853214)], [(38853353, 38853472)], [(38855530, 38855611)], [(38855700, 38855757)], [(38857795, 38857952)], [(38858148, 38858212)], [(38858320, 38858416)], [(38858687, 38858777)], [(38860611, 38860612)]], ['VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARTRUE']))
    """
    if splitters is None:
        return [[x] for x in ranges], varflags
    results=[x._vals for x in IntervalSet(ranges).split(splitters)]
    vf=[]
    for i, iv in enumerate(results):
        for j, r in enumerate(ranges):
            if overlaps(iv[0][0], iv[0][1], r[0], r[1]):
                vf.append(varflags[j])
                break
    assert len(results) == len(vf)
    return results, vf

def get_ranges(last, vstart, vend, exon_starts, exon_ends, chrom=1): # NOTE: new model version
    """
    >>> get_ranges(61018, 62029, 62029, (
    ... 60174, 60370, 60665, 60925, 62029, 62216, 62453,
    ... 62675, 63052, 63398, 63652, 63868, 64512, 64764,
    ... 65018, 65671), (60281,
    ... 60565, 60808, 61033, 62134, 62379, 62587, 62824,
    ... 63209, 63559, 63779, 64102, 64691, 64946, 65084,
    ... 65985))
    ([(61018, 61033)], 61018, ['VARFALSE'])

    >>> get_ranges(617350, 617346, 617369, (617350, 617350), (617400, 617400))
    ([], 617369, [])
    
    >>> get_ranges(0, 1, 1, (0, 20), (10, 30)) # needs varflag
    ([(0, 1)], 0, ['VARTRUE'])

    >>> get_ranges(0, 3, 3, (0, 20), (10, 30)) # second region needs varflag
    ([(0, 2), (2, 3)], 0, ['VARFALSE', 'VARTRUE'])

    >>> get_ranges(3, 6, 6, (0, 20), (10, 30)) # second region needs varflag
    ([(3, 5), (5, 6)], 3, ['VARFALSE', 'VARTRUE'])

    >>> get_ranges(6, 7, 7, (0, 20), (10, 30)) # needs varflag
    ([(6, 7)], 6, ['VARTRUE'])

    >>> get_ranges(7, 9, 9, (0, 20), (10, 30)) # needs varflag for second region only
    ([(7, 8), (8, 9)], 7, ['VARFALSE', 'VARTRUE'])

    >>> get_ranges(38865400, 38865404, 38865425, (38865242, 38865242), (38865601, 38865601))
    ([(38865400, 38865403), (38865403, 38865425)], 38865425, ['VARFALSE', 'VARTRUE'])
    
    >>> get_ranges(38865405, 38865404, 38865425, (38865242, 38865242), (38865601, 38865601))
    ([], 38865425, [])

    >>> get_ranges(61018, 61990, 62001, (60925, 62000), (61033, 62040))
    ([(61018, 61033), (62000, 62001)], 62001, ['VARFALSE', 'VARTRUE'])

    >>> get_ranges(61018, 62023, 62030, (60925, 62000), (61033, 62040)) # this situation in the future could possibly become (62022, 62030) with the vartrue flag
    ([(61018, 61033), (62000, 62022), (62022, 62030)], 62030, ['VARFALSE', 'VARFALSE', 'VARTRUE'])
    
    >>> get_ranges(62000, 62023, 62050, (60925, 62000, 62045), (61033, 62040, 62060)) # this situation in the future could possibly become (62022, 62030) with the vartrue flag
    ([(62000, 62022), (62022, 62040), (62045, 62050)], 62050, ['VARFALSE', 'VARTRUE', 'VARTRUE'])

    >>> get_ranges(62000, 62023, 62070, (60925, 62000, 62045), (61033, 62040, 62060)) # this situation in the future could possibly become (62022, 62030) with the vartrue flag
    ([(62000, 62022), (62022, 62040), (62045, 62060)], 62070, ['VARFALSE', 'VARTRUE', 'VARTRUE'])

    >>> get_ranges(56, 95, 95, range(0, 1000, 10), range(5, 1000, 10))
    ([(60, 65), (70, 75), (80, 85), (90, 94), (94, 95)], 60, ['VARFALSE', 'VARFALSE', 'VARFALSE', 'VARFALSE', 'VARTRUE'])

    >>> get_ranges(0, 10, 10, range(0, 100, 10), range(5, 100, 10))
    ([(0, 5)], 0, ['VARFALSE'])

    """

    f4=open('deletioncut.txt','a') #code removed by deletions

    varflag=[]
 
    assert last >= exon_starts[0]
    assert vstart <= exon_ends[-1]
    assert all(s < e for s, e in zip(exon_starts, exon_ends))

    istart = bisect_left(exon_starts, last) - 1
    if istart == -1: istart = 0
    [(0, 6), (10, 11)]

    assert exon_starts[istart] <= last, (exon_starts[istart], last, istart)
    if exon_ends[istart] <= last:
        istart += 1
        if istart < len(exon_starts): # in case last and exon ends are equal, it will loop through again, but I want that last loop
            last = exon_starts[istart]
    start = last

    if vstart < vend: # moved here because there are variants in UTRs that do not exist in coding exon space
        last = vend
    ranges = []
    while start < vstart and istart < len(exon_starts): #<= lets it capture 0 length regions, so I removed it and the +1 allows it to make 1 bp regions when two variants are right next to one another
        ranges.append((start, exon_ends[istart])) #removed +1 from exon_ends[istart] + 1, because IntervalSet is already in 0-based half-open format
        istart += 1
        varflag.append("VARFALSE") # unless using vstart, there is no variant
        try: 
            if exon_starts[istart] > vstart and exon_starts[istart] < vend and ranges[-1][1] < vstart:
                ranges.append((exon_starts[istart], vend))
                varflag.append("VARTRUE")
                break
        except IndexError:
            pass
        if ranges[-1][1] >= vstart: # equal to is now possible, since we are including variant start+1 and ranges are in 0-based half-open
            if ranges[-1][0]-(vstart-1)==0:
                ranges[-1] = (ranges[-1][0], vstart)
                varflag[-1]="VARTRUE"
                break
            ranges[-1] = (ranges[-1][0], vstart-1) #removed +1 from vstart + 1, because IntervalSet is already in 0-based half-open format
            varflag.append("VARTRUE") #variant contained at end coordinate = TRUE; this indicates the region contains the variant, therefore should be considered 0 bp, get a 0 coverage and a 0 cpg score
            if exon_ends[istart-1] < vend:
                varflag.append("VARTRUE")
                if exon_ends[-1] < vend:
                    ranges.append((vstart-1, exon_ends[istart-1])); ranges.append((exon_starts[istart], exon_ends[-1]))
                    break
                ranges.append((vstart-1, exon_ends[istart-1])); ranges.append((exon_starts[istart], vend))
                break
            ranges.append((vstart-1, vend))
            break
        start = exon_starts[istart]

    return ranges, last, varflag

def get_ranges_w_variant(last, vstart, vend, exon_starts, exon_ends, chrom=1): # NOTE: the previous version of the model where end coordinate contains variant
    """
    >>> get_ranges_w_variant(61018, 62029, 62029, (
    ... 60174, 60370, 60665, 60925, 62029, 62216, 62453,
    ... 62675, 63052, 63398, 63652, 63868, 64512, 64764,
    ... 65018, 65671), (60281,
    ... 60565, 60808, 61033, 62134, 62379, 62587, 62824,
    ... 63209, 63559, 63779, 64102, 64691, 64946, 65084,
    ... 65985))
    ([(61018, 61033)], 61018, ['VARFALSE'])

    >>> get_ranges_w_variant(38865400, 38865404, 38865425, (38865242, 38865242), (38865601, 38865601))
    ([(38865400, 38865404)], 38865425, ['VARTRUE'])
    
    >>> get_ranges_w_variant(38865405, 38865404, 38865425, (38865242, 38865242), (38865601, 38865601))
    ([], 38865425, [])

    >>> get_ranges_w_variant(617350, 617346, 617369, (617350, 617350), (617400, 617400))
    ([], 617369, [])

    >>> get_ranges_w_variant(61018, 61990, 62001, (60925, 62000), (61033, 62040))
    ([(61018, 61033)], 62001, ['VARFALSE'])

    >>> get_ranges_w_variant(61018, 62023, 62030, (60925, 62000), (61033, 62040)) #varflag for last one only
    ([(61018, 61033), (62000, 62023)], 62030, ['VARFALSE', 'VARTRUE'])

    >>> get_ranges_w_variant(56, 95, 95, range(0, 1000, 10), range(5, 1000, 10)) #varflag for last one only
    ([(60, 65), (70, 75), (80, 85), (90, 95)], 60, ['VARFALSE', 'VARFALSE', 'VARFALSE', 'VARTRUE'])

    >>> get_ranges_w_variant(1, 10, 10, range(0, 100, 10), range(5, 100, 10))
    ([(1, 5)], 1, ['VARFALSE'])

    >>> get_ranges_w_variant(0, 10, 10, range(0, 100, 10), range(5, 100, 10))
    ([(0, 5)], 0, ['VARFALSE'])

    >>> get_ranges_w_variant(50, 59, 59, (50, 61), (60, 70))
    ([(50, 59)], 50, ['VARTRUE'])

    >>> get_ranges_w_variant(1562576, 1562675, 1562675,
    ... (1560665, 1560925, 1562029, 1562216, 1562453, 1562675, 1563052, 1563398, 1563652, 1563868, 1564512, 1564764, 1565018, 1565671),
    ... (1560808, 1561033, 1562134, 1562379, 1562587, 1562824, 1563209, 1563559, 1563779, 1564102, 1564691, 1564946, 1565084, 1565985))
    ([(1562576, 1562587)], 1562576, ['VARFALSE'])
    """

    f4=open('deletioncut.txt','a') #code removed by deletions

    varflag=[]
 
    assert last >= exon_starts[0]
    assert vstart <= exon_ends[-1]
    #assert vstart >= last, (vstart, last, exon_starts) # deletion can overlap UTR/intron, but this is controlled for in exac-regions.py
    assert all(s < e for s, e in zip(exon_starts, exon_ends))

    istart = bisect_left(exon_starts, last) - 1
    if istart == -1: istart = 0
    [(0, 6), (10, 11)]

    assert exon_starts[istart] <= last, (exon_starts[istart], last, istart)
    if exon_ends[istart] <= last:
        istart += 1
        if istart < len(exon_starts): # in case last and exon ends are equal, it will loop through again, but I want that last loop
            last = exon_starts[istart]
    start = last

    if vstart < vend: # moved here because there are variants in UTRs that do not exist in coding exon space
        last = vend
    ranges = []; writedel=True
    while start < vstart and istart < len(exon_starts): #<= lets it capture 0 length regions, so I removed it and the +1 allows it to make 1 bp regions when two variants are right next to one another
        ranges.append((start, exon_ends[istart])) #removed +1 from exon_ends[istart] + 1, because IntervalSet is already in 0-based half-open format
        if writedel and vstart < vend and vend > exon_starts[istart]:
            f4.write("\t".join(map(str,[chrom,vstart,last]))+"\n") # removed by deletion, but only if it is within an exon
            writedel=False
        istart += 1
        varflag.append("VARFALSE") # unless using vstart, there is no variant
        if ranges[-1][1] >= vstart: # equal to is now possible, since we are including variant start+1 and ranges are in 0-based half-open
            ranges[-1] = (ranges[-1][0], vstart) #removed +1 from vstart + 1, because IntervalSet is already in 0-based half-open format
            varflag[-1]="VARTRUE" #variant contained at end coordinate = TRUE
            break
        start = exon_starts[istart]

    return ranges, last, varflag

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
    p = subprocess.Popen("tabix {path}/Panel.chr{chrom}.coverage.txt.gz {chrom} | cut -f 2,{coli} ".format(**locals()), # {path}/exacv2.chr{chrom}.cov.txt.gz
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
    f1 = open("selfchaincut.txt","a")
    f2 = open("segdupscut.txt","a")
    f3 = open("coveragecut.txt","a")
    for a in args:
        assert a[0] == "|", ("expecting a tabix query so we can handle chroms correctly", a)
    
        # any file that gets sent in will be used to split regions (just like
        # low-coverage). For example, we split on self-chains as well.
#TODO: comment this block if you don't want any filtering by self-chains or segdups
        for toks in (x.strip().split("\t") for x in ts.nopen(a)): # adds self chains and segdups to splitters list, so that exons can be split, and they are removed from CCRs
            s, e = int(toks[1]), int(toks[2])
            split_iv.add((s, e))
            if len(toks) > 3:
                f1.write("\t".join(toks)+"\n") # self chain
            else:
                f2.write("\t".join(toks)+"\n") # segdups
                

    for toks in (x.rstrip('\r\n').split("\t") for x in ts.nopen(gtf) if x[0] != "#"):
        if toks[2] not in("CDS", "stop_codon") or toks[1] not in("protein_coding"): continue
        #if toks[0] != "1": break
        start, end = map(int, toks[3:5])
        gene = toks[8].split('gene_name "')[1].split('"', 1)[0]
        assert start <= end, toks
        key = toks[0], gene

        cutoff = 0.3

        # find sections of exon under certain coverage.
#TODO: comment this if we don't want coverage cutoff filtering
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

        genes[key].add([(start-1, end)]) # converts GTF exon coordinates to BED format (subtracts 1 from exon start)
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
    return any(c in csq['Consequence'] for c in ('stop_gained', 'stop_lost', 'start_lost', 'initiator_codon_variant', 'rare_amino_acid_variant', 'missense_variant', 'protein_altering_variant', 'frameshift_variant', 'inframe_insertion', 'inframe_deletion')) \
    or (('splice_donor_variant' in csq['Consequence'] or 'splice_acceptor_variant' in csq['Consequence'] or '5_prime_UTR_variant' in csq['Consequence'] or '3_prime_UTR_variant' in csq['Consequence']) and 'coding_sequence_variant' in csq['Consequence'])

def cg_content(seq):
    if len(seq) == 0: return 0.0
    return 2.0 * seq.count('CG') / len(seq)
