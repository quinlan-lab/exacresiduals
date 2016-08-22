import os
import sys
from collections import defaultdict, OrderedDict
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
import subprocess
import numpy as np
import toolshed as ts
import seaborn as sns
sns.set_style('whitegrid')

from cyvcf2 import VCF
import utils as u
import pyinter
from bw import BigWig

region = sys.argv[1]
if len(sys.argv) > 2:
    assert sys.argv[2].isdigit() and sys.argv[3].isdigit()
    region = "%s:%s-%s" % (sys.argv[1], sys.argv[2], sys.argv[3])
elif "\t" in region:
    toks = region.split("\t")
    region = "%s:%s-%s" % tuple(toks[:3])

def read_variants(region, path="data/ExAC.r0.3.sites.vep.vcf.gz"):
    """
    read ExAC coverage from a single chrom into a numpy array. If no length is
    given, just use the one length from chrom 1.
    path is expected to contain Panel.chr*
    info field is the column to pull
    """
    
    # just extract the position (2) and the requested column
    vcf = VCF(path)
    
    filters=[]
    a=vcf.raw_header # the following block of code gets the filter list
    b=a.split("\n")
    for i in b:
        if '##FILTER' in i:
            filters.append(i.split('ID=')[1].split(',')[0])
 
    var = defaultdict(list)
    j = 0
    kcsq = vcf["CSQ"]["Description"].split(":")[1].strip(' "').split("|")
    for v in vcf(region):
        csqs = [dict(zip(kcsq, c.split("|"))) for c in v.INFO['CSQ'].split(",")]
        for csq in (c for c in csqs if c['BIOTYPE'] == 'protein_coding'):
            if csq['Feature'] == '' or csq['EXON'] == '' or csq['cDNA_position'] == '': continue
            for c in csq['Consequence'].split('&'):
                if c in ('synonymous_variant'):
                    var['syn'].append(v.start)
                    break
        if v.FILTER is not None and any(u.isfunctional(csq) for csq in csqs):
            var[v.FILTER].append(v.start)
        j += 1
        #if j > 100000: break
    assert j > 0, ("no values found for", chrom, path)
    return var, filters

def read_gerp(region, sends, path='/scratch/ucgd/lustre/u1021864/serial/hg19.gerp.bw'):
    gerp = BigWig(path)
    exongerp=[]
    chrom, se = region.split(":")
    s, e = map(int, se.split("-"))
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    for key in sends:
        for (exs, exe) in zip(sends[key][0], sends[key][1]):
             exongerp.extend(np.frombuffer(gerp.values(chrom, int(exs)-1, int(exe)), dtype='f'))

    return np.frombuffer(gerp.values(chrom, int(s)-1, int(e)), dtype='f'), exongerp

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
    totlen = 0
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
    ints={}
    lens=pyinter.IntervalSet()
    for tr, ivset in transcripts.iteritems():
        sends = sorted(list(ivset))
        iset=pyinter.IntervalSet(pyinter.closedopen(x.lower_value,x.upper_value) for x in sends)
        lens = lens.union(iset)
        ss, es = [x.lower_value for x in sends], [x.upper_value for x in sends]
        ints[tr] = (ss,es)
    totlen = sum(x.upper_value-x.lower_value for x in lens)
    return ints, set(names), set(ids), set(trs), totlen

def read_pfam(path):
    tracks = defaultdict(pyinter.IntervalSet)
    pids, trs, ids = [], [], []
    for toks in (x.rstrip('\r\n').split() for x in ts.nopen(path) if x[0] != "#"):
        start, end = map(int, toks[1:3])
        assert start <= end, toks
        pid = toks[10].split(';',1)[0].strip('"') #pfamA_id
        tracks[pid].add(pyinter.closedopen(start, end))
        ids.append(toks[12].split(';',1)[0].strip('"')) #gene_name
        trs.append(toks[14].split(';',1)[0].strip('"')) #transcript_id
    ints={}
    for pid, ivset in tracks.iteritems():
        sends = sorted(list(ivset))
        ss, es = [x.lower_value for x in sends], [x.upper_value for x in sends]
        ints[pid] = (ss,es)

    return ints

def read_repeats(path,keyname):
    tracks = defaultdict(pyinter.IntervalSet)
    for toks in (x.rstrip('\r\n').split() for x in ts.nopen(path) if x[0] != "#"):
        start, end = map(int, toks[1:3])
        assert start <= end, toks
        tracks[keyname].add(pyinter.closedopen(start, end))
    ints={}
    for pid, ivset in tracks.iteritems():
        sends = sorted(list(ivset))
        ss, es = [x.lower_value for x in sends], [x.upper_value for x in sends]
        ints[pid] = (ss,es)

    return ints

#sort -k1,1 -k2,2n $DATA/hgsegmental.bed | uniq | bedtools merge | bgzip -c > data/segmental.bed.gz

gd=OrderedDict()
gs,ge={},{}
keys=[]
f, axarr = plt.subplots(3, sharex=True)

sends, names, ids, trs, totlen = read_exons("| tabix /scratch/ucgd/lustre/u1021864/serial/Homo_sapiens.GRCh37.75.gtf.gz {region}".format(region=region))
gd.update(sends)
keys.extend(sends.keys())

s, e, cov = read_coverage(region)
axarr[0].plot(range(s, e + 1), cov)
ymin,ymax=axarr[0].get_ylim()[0]-.05,axarr[0].get_ylim()[1]+.05
axarr[0].set_ylim(ymin,ymax)

gerp, exongerp = read_gerp(region, sends)
axarr[1].plot(range(s, e + 1), gerp)
ymin,ymax=-12.36,6.18 
#ymin,ymax=axarr[1].get_ylim()[0]-.05,axarr[1].get_ylim()[1]+.05axarr[1].set_ylim(ymin,ymax)

sends = read_pfam("| tabix data/pfam.bed.gz {region}".format(region=region))
gd.update(sends)
keys.extend(sends.keys())

sends = read_repeats("| tabix data/self-chains.gt90.bed.gz {region}".format(region=region),'selfchain')
gd.update(sends)
keys.extend(sends.keys())

sends = read_repeats("| tabix data/segmental.bed.gz {region}".format(region=region), 'segdup')
gd.update(sends)
keys.extend(sends.keys())

var, filters = read_variants(region)
var2 = OrderedDict(sorted(var.items(), key=lambda t: t[0]))

"""
calculating VQSR Tranche densities
"""

vqsr={}
for i in var2:
    if i.startswith('VQSR'):
        vqsr[i]=int(1/(len(var2[i])/float(totlen)))
vqsrlen=len(vqsr)
gd.update(var2)
keys.extend(var2.keys())
markers = ['bo','ro','go','yo','mo','co','ko']
j = 0

axarr[0].set_title("coverage plot -- sum(exonic cov): %.1f" % (cov.sum()))
axarr[1].set_title("gerp plot -- mean(exonic gerp): %.1f" % (np.mean(exongerp)))
densities=[]
for i, k in vqsr.items():
    densities.append(i); densities.append(k)
strcat="; ".join(['%s density: 1/%s' for i in vqsr.keys()]) % tuple(densities)
strcat="%s/%s %s -- syn density: 1/%i\n" % ("|".join(names), "|".join(ids), region, int(1/(len(var['syn'])/float(totlen)))) + strcat
axarr[2].set_title(strcat.replace('VQSRTranche',''))
plt.subplots_adjust(hspace=0.5)
#bbox = f.get_tightbbox(f.canvas.get_renderer())

for ind, key in enumerate(gd):
    if key.startswith('VQSR'):
        marker = 's'
        j+=1
        color = (0,1/float(vqsrlen)*j,0)
    elif key == 'syn':
        marker = 'o'
        color = 'c'
    elif key in filters and not key.startswith('VQSR'): 
        marker = 'o'
        color = 'k'
    else:
        marker = ''
        color = 'k'
        if key.startswith('ENST'):
            color = 'b'
    if marker != '':
        axarr[2].plot(gd[key], np.zeros(len(gd[key])) + (ind+1)/1., marker=marker, color=color, label = key, ls='none')
    else:
        for k, (exs, exe) in enumerate(zip(gd[key][0], gd[key][1])):
            axarr[2].plot([exs, exe], [ind+1/1., ind+1/1.], ls='-', color=color, lw=3)

axarr[2].set_yticks(np.arange(1,(ind+2)/1))
fp=matplotlib.font_manager.FontProperties(family='sans-serif', style='normal', variant='normal', weight='normal', stretch='normal', size='xx-small')
axarr[2].set_yticklabels(keys, fontproperties=fp)
ymin,ymax=axarr[2].get_ylim()[0]-.25,axarr[2].get_ylim()[1]+.25
axarr[2].set_ylim(ymin,ymax)

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
plt.savefig('figs/' + region.replace(":", "-") + ".pdf", format='pdf', bbox_inches='tight')
plt.close()
