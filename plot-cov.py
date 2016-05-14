import os
import sys
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
import subprocess
import numpy as np
import toolshed as ts
import seaborn as sns
sns.set_style('whitegrid')


region = sys.argv[1]
if len(sys.argv) > 2:
    assert sys.argv[2].isdigit() and sys.argv[3].isdigit()
    region = "%s:%s-%s" % (sys.argv[1], sys.argv[2], sys.argv[3])
elif "\t" in region:
    toks = region.split("\t")
    print(toks)
    region = "%s:%s-%s" % tuple(toks[:3])

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

    names = []
    trs, ids = [], []
    for toks in (x.rstrip('\r\n').split("\t") for x in ts.nopen(gtf) if x[0] != "#"):
        if toks[2] not in("CDS", "stop_codon") or toks[1] not in("protein_coding"): continue
        #if toks[0] != "1": break
        start, end = map(int, toks[3:5])
        assert start <= end, toks
        transcript = toks[8].split('transcript_id "')[1].split('"', 1)[0]
        transcripts[transcript].add(pyinter.closedopen(start-1, end))

        names.append(toks[8].split('transcript_name "')[1].split('"', 1)[0].rsplit("-", 1)[0])
        ids.append(toks[8].split('gene_id "')[1].split('"', 1)[0])
        trs.append(toks[8].split('transcript_id "')[1].split('"', 1)[0])

    # sort by start so we can do binary search.
    # TODO: need to remove overlapping exons so we don't double-count
    transcripts = dict((k, sorted(v)) for k, v in transcripts.iteritems())
    #ends = dict((k, sorted(v)) for k, v in ends.iteritems())
    starts, ends = {}, {}
    for tr, ivset in transcripts.iteritems():
        sends = sorted(list(ivset))
        starts[tr] = [x.lower_value for x in sends]
        ends[tr] = [x.upper_value for x in sends]
    return starts, ends, set(names), set(ids), set(trs)

import pyinter


s, e, cov = read_coverage(region)
print(len(range(s, e)), len(cov))
plt.plot(range(s, e + 1), cov)
starts, ends, names, ids, trs = read_exons("| tabix /scratch/ucgd/lustre/u1021864/serial/Homo_sapiens.GRCh37.75.gtf.gz {region}".format(region=region))


for i, tr in enumerate(starts, start=1):
    for k, (exs, exe) in enumerate(zip(starts[tr], ends[tr])):
        plt.plot([exs, exe], [-0.08 * i, -0.08 * i], 'k-', lw=3)
        if k == 0:
            plt.text(exe + 50, -0.08 * i + 0.01, tr)

plt.title("%s/%s %s -- sum(cov): %.1f" % ("|".join(names), "|".join(ids),
    region, cov.sum()))
plt.xlim(s, e)
plt.draw()
ticks, labels = plt.xticks()
if len(ticks) > 4:
    ticks = ticks[::2]
    labels = labels[::2]
plt.xticks(ticks, labels)

from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
ax = plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

plt.xlim(s, e)

plt.savefig('figs/' + region.replace(":", "-") + ".png")
plt.close()
