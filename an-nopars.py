from cyvcf2 import VCF
import tabix

#pars hard to find, but at: http://genome-test.soe.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=map&hgta_track=par&hgta_table=par (agrees with NCBI gtf, I checked)

vcf = VCF('data/gnomad-vep-vt.vcf.gz')
an=0
tb = tabix.open('par.bed.gz')
for v in vcf:
    if v.CHROM == "X":
        r = tb.querys(v.CHROM+":"+str(v.POS-1)+"-"+str(v.POS))
        l = [i for i in r]
        if not l:
            an = max(v.INFO['AN'], an)
print "Allele Number Max"
print an
