from cyvcf2 import VCF
import argparse
import numpy as np
import subprocess as sp
import utils as u


#intersect variants, add exacvars/gnomad vars to .sh script, then just tabix query gnomad for variants with greatest change and run statistics.

parser=argparse.ArgumentParser()
parser.add_argument("-e", "--exac", help="exac variants")
parser.add_argument("-g", "--gnomad", help="gnomad variants")
parser.set_defaults(exac = '/uufs/chpc.utah.edu/common/home/u1021864/analysis/exacresiduals/data/ExAC.r1.vt.vep.vcf.gz')
parser.set_defaults(gnomad = '/uufs/chpc.utah.edu/common/home/u1021864/analysis/exacresiduals/data/gnomad.exomes.r2.0.1.sites.vep.vt.vcf.gz')
args=parser.parse_args()
exac=args.exac
gnomad=args.gnomad

exac = VCF(exac)

gnomad = VCF(gnomad)
kcsq = gnomad["CSQ"]["Description"].split(":")[1].strip(' "').split("|")

f=open("diff.txt","r")
f2=open("results.txt", "w")

prevpos=-1; idx=0; varct, sct = 0, 0; prevvar=""
for line in f:
    fields = line.strip().split('\t')
    chrom=fields[0]; start=fields[1]; end=fields[2]
    position=chrom+":"+str(int(start)+1)+"-"+end
    for gnovar in gnomad(position):
        if str(gnovar)==str(prevvar):
            continue
        prevvar=gnovar
        info = gnovar.INFO
        try:
            csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]
        except KeyError:
            continue
        if prevpos == gnovar.POS:
            idx+=1
        else:
            idx=0
            prevpos = gnovar.POS
        for csq in (c for c in csqs if c['BIOTYPE'] == 'protein_coding'):
            if u.isfunctional(csq):
                varct+=1
                #f2.write(str(gnovar))
                ac = info['AC']
                if not isinstance(info['AC'], (int, long)):
                    #print gnovar, idx
                    ac = ac[idx] 
                if ac == 1:
                    sct+=1
                    #f2.write(str(gnovar))
                break
            
f2.write(str(sct/float(varct)) + "\n")

def query(regions,hotspot,wo):
    p1 = sp.Popen('zcat '+regions, shell = True, stdout = sp.PIPE)
    l=['bedtools', 'intersect', '-a', 'stdin', '-b', hotspot]
    if wo:
        l.append('-wo')
    p2 = sp.Popen(l, stdin = p1.stdout, stdout = sp.PIPE)
    output,error = p2.communicate()
    p1.kill()
    return output.strip()

#result = query(exac, gnomad, 1)
#f.write(result)
f.close(); f2.close()
