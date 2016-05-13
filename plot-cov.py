import os
import sys
from collections import defaultdict

from matplotlib import pyplot as plt
import subprocess
import numpy as np
import toolshed as ts
import seaborn as sns


region = sys.argv[1]



def read_coverage(region, cov=10, path="~u6000771/Data/ExAC-coverage/"):
    """
    read ExAC coverage from a single chrom into a numpy array. If no length is
    given, just use the one length from chrom 1.
    path is expected to contain Panel.chr*
    cov is the column to pull
    """

    cols = "chrom	pos	mean	median	1	5	10	15	20	25	30	50 100".split()
    coli = cols.index(str(cov)) + 1


    chrom, se = region.split(":")
    s, e = map(int, se.split("-"))
    length = e - s + 1

    # just extract the position (2) and the requested column
    p = subprocess.Popen("tabix {path}/Panel.chr{chrom}.coverage.txt.gz {region} | cut -f 2,{coli} ".format(**locals()),
            stdout=subprocess.PIPE, stderr=sys.stderr,
            shell=True,
            executable=os.environ.get("SHELL"))

    cov = np.zeros(length, dtype=np.float32)
    j = 0
    for line in p.stdout:
        pos, val = line.split()
        cov[int(pos)-s] = float(val)
        j += 1
        #if j > 100000: break
    assert j > 0, ("no values found for", chrom, path)
    p.wait()
    if p.returncode != 0:
        raise Exception("bad: %d", p.returncode)
    return s, e, cov

def read_exons(gtf):
    transcripts = defaultdict(pyinter.IntervalSet)

    for toks in (x.rstrip('\r\n').split("\t") for x in ts.nopen(gtf) if x[0] != "#"):
        if toks[2] not in("CDS", "stop_codon") or toks[1] not in("protein_coding"): continue
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
    fstarts, fends = [], []
    for tr in starts:
        fstarts.extend(starts[tr])
        fends.extend(ends[tr])
    return fstarts, fends

import pyinter


s, e, cov = read_coverage(region)
print(len(range(s, e)), len(cov))
plt.plot(range(s, e + 1), cov)
starts, ends = read_exons("| tabix /scratch/ucgd/lustre/u1021864/serial/Homo_sapiens.GRCh37.75.gtf.gz {region}".format(region=region))

for exs, exe in zip(starts, ends):
    plt.plot([exs, exe], [-0.1, -0.1], 'k-', lw=4)

plt.xlim(s, e)
plt.title("%s -- sum(cov): %.1f" % (region, cov.sum()))
plt.show()
