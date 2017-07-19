import subprocess

f1=open('flatexome.bed','r')
f2=open('cutregions.txt','r')
f3=open('results/2016_12_10/resids.txt')
prevgene=None
gbed=''
for line in f:
    fields=line.strip().split('\t')
    chrom=fields[0];start=fields[1];end=fields[2];transcript=fields[3];gene=fields[4]
    if gene==prevgene or prevgene is None:
        gbed+="\t".join([chrom,start,end])+"\n"
    else:
        p=subprocess.Popen(['bedtools merge'],shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
        output,err = p.communicate(gbed)
        gbed="\t".join([chrom,start,end])+"\n"
        if prevgene:
            for line in output.strip().split("\n"):
                print line + "\t" + prevgene
        else:
            for line in output.strip().split("\n"):
                print line + "\t" + gene
    prevgene=gene
