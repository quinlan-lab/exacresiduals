import sys
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-b", "--basepair", help="use this argument if you want a match by basepair", action="store_true", default=False)
parser.add_argument("-g", "--genes", help="use this argument if you want a match by number of genes", action="store_true", default=False)
args=parser.parse_args()

if args.basepair:
    f = open('/uufs/chpc.utah.edu/common/home/u1021864/analysis/regions/topresid.txt','r')
    totbp = 0
    for line in f:
        fields = line.strip().split("\t")
        totbp += int(fields[2]) - int(fields[1])
    f.close()
    f = open('/uufs/chpc.utah.edu/common/home/u1021864/analysis/exacresiduals/results/2016_06_15/exonicresiduals.txt', 'r')
    list_ = []
    f.readline()
    for line in f:
        fields = line.strip().split("\t")
        list_.append(fields)
    bp = 0 
    list_ = sorted(list_,key = lambda x: abs(float(x[10])))
    for line in list_:
        if bp > totbp:
            break
        bp += int(line[2]) - int(line[1])
        print "\t".join(line)

if args.genes:
    f = open('/uufs/chpc.utah.edu/common/home/u1021864/analysis/regions/topresid.txt','r')
    topgeneset=set()
    for line in f:
        fields = line.strip().split("\t")
        topgeneset.add(fields[3])
    f.close()
    f = open('/uufs/chpc.utah.edu/common/home/u1021864/analysis/exacresiduals/results/2016_06_15/exonicresiduals.txt', 'r')
    list_ = []
    f.readline()
    for line in f:
        fields = line.strip().split("\t")
        list_.append(fields)
    totgenes=len(topgeneset)
    list_ = sorted(list_,key = lambda x: abs(float(x[10])))
    geneset=set()
    for line in list_:
        geneset.add(line[3])
        genes=len(geneset)
        if genes > totgenes:
            break
        print "\t".join(line)
