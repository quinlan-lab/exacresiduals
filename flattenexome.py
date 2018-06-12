import subprocess
import sys
from itertools import groupby
from operator import itemgetter

f=open(sys.argv[1],'r') # codingtranscriptome.bed
genes = []
for line in f:
    fields=line.strip().split('\t')
    chrom=fields[0];start=int(fields[1]);end=int(fields[2]);transcript=fields[3];gene=fields[4]
    genes.append((chrom, start, end, transcript, gene))
sorter = itemgetter(4,0,1,2)
grouper = itemgetter(4)
gbed=''
for key, grp in groupby(sorted(genes, key = sorter), grouper):
    grp=list(grp)
    gene=grp[0][-1]
    for i, elem in enumerate(grp):
        chrom=grp[i][0]; start=str(grp[i][1]); end=str(grp[i][2])
        gbed+="\t".join([chrom,start,end])+"\n"
    p=subprocess.Popen(['bedtools merge'],shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    output,err = p.communicate(gbed)
    for line in output.strip().split("\n"):
        print line + "\t" + gene
    gbed=''
