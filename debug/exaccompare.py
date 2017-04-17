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

f1=open("pathodiff.txt", "r")
f2=open("benigndiff.txt", "r")
f3=open("pathoresults.txt", "w")
f4=open("benignresults.txt", "w")

def stats(infile, outfile):
    prevpos = -1; idx = 0; varct, sct, passct, uct, dct = 0.0, 0.0, 0.0, 0.0, 0.0; prevvar = None; prevline = None
    upct,dpct=[],[]
    for line in infile:
        fields = line.strip().split('\t')
        chrom=fields[0]; start=fields[1]; end=fields[2]
        position=chrom+":"+str(int(start)+1)+"-"+end
        if prevline:
            if line == prevline:
                continue
        prevline=line
        for gnovar in gnomad(position):
            if prevpos == gnovar.POS:
                idx+=1
            else:
                idx=0
                prevpos = gnovar.POS
            #print line, gnovar, idx
            if prevvar:
                newkey=gnovar.CHROM+str(gnovar.POS)+str(gnovar.ID)+str(gnovar.REF)+str(gnovar.ALT)
                oldkey=prevvar.CHROM+str(prevvar.POS)+str(prevvar.ID)+str(prevvar.REF)+str(prevvar.ALT)
                if newkey==oldkey:
                    continue
            prevvar=gnovar
            info = gnovar.INFO
            try:
                csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]
            except KeyError:
                continue
            for csq in (c for c in csqs if c['BIOTYPE'] == 'protein_coding'):
                if u.isfunctional(csq):
                    for exacvar in exac(position):
                        if str(gnovar.CHROM)+str(gnovar.POS)+str(gnovar.REF)+str(gnovar.ALT)==str(exacvar.CHROM)+str(exacvar.POS)+str(exacvar.REF)+str(exacvar.ALT):
                            if gnovar.FILTER!=exacvar.FILTER:
                                if (gnovar.FILTER is None or gnovar.FILTER in ["PASS", "SEGDUP"]) and not (exacvar.FILTER is None or exacvar.FILTER in ["PASS", "SEGDUP"]):
                                    upct.append(float(fields[4])-float(fields[9]))
                                    #upct+=float(fields[4])
                                    uct+=1 # quality upgrade to PASS
                                    # outfile.write(str(gnovar))
                                    # outfile.write(str(exacvar)+"\n")
                                if not (gnovar.FILTER is None or gnovar.FILTER in ["PASS", "SEGDUP"]) and (exacvar.FILTER is None or exacvar.FILTER in ["PASS", "SEGDUP"]):
                                    dpct.append(float(fields[4])-float(fields[9]))
                                    #dpct+=float(fields[4])
                                    dct+=1 # quality downgrade from PASS
                                    # outfile.write(str(gnovar))
                                    # outfile.write(str(exacvar+"\n"))
                    varct+=1 #func var
                    #outfile.write(str(gnovar))
                    ac = info['AC']
                    if not isinstance(info['AC'], (int, long)):
                        ac = ac[idx] 
                    if ac == 1:
                        sct+=1 # singleton
                        #outfile.write(str(gnovar))
                    if not (gnovar.FILTER is None or gnovar.FILTER in ["PASS", "SEGDUP"]):
                        passct+=1 #non-PASS gnovar
                    break

    outfile.write("singleton: " + "%.3f" % (100*sct/float(varct)) + "\n")
    outfile.write("non-pass gnomad: " + "%.3f" % (100*passct/float(varct)) + "\n")
    outfile.write("downgrade filter from exac to gnomad: " + "%.3f" % (100*dct/float(varct)) + "\n")
    outfile.write("upgrade to pass filter from exac to gnomad: " + "%.3f" % (100*uct/float(varct)) + "\n")
    outfile.write("average percentile difference of downgraded variants: " + "%.3f" % (sum(dpct)/float(1 if dct == 0 else dct)) +  "\n")
    outfile.write("average percentile difference of upgraded variants: " + "%.3f" % (sum(upct)/float(1 if uct == 0 else uct)) +  "\n")
    #outfile.write("percentiles of upgraded variants: " + str(upct) + "\n")
    #outfile.write("percentiles of downgraded variants: " + str(dpct) + "\n")

stats(f1,f3)
stats(f2,f4)

f1.close(); f2.close(); f3.close(); f4.close()

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
