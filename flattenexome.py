import subprocess

f=open('codingtranscriptome.bed','r')
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
p=subprocess.Popen(['bedtools merge'],shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
output,err = p.communicate(gbed)
for line in output.strip().split("\n"):
    print line + "\t" + prevgene

