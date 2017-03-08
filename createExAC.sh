# get ExAC v2.0 VCF, coverage files
wget -P $DATA https://storage.googleapis.com/gnomad-public/release-170228/vcf/exomes/gnomad.exomes.r2.0.1.sites.vcf.gz
for i in {1..22} X Y; do wget -P $HOME/analysis/exacresiduals/data https://storage.googleapis.com/gnomad-public/release-170228/coverage/exomes/exacv2.chr$i.cov.txt.gz; done
for i in {1..22} X Y; do wget -P $HOME/analysis/exacresiduals/data https://storage.googleapis.com/gnomad-public/release-170228/coverage/exomes/exacv2.chr$i.cov.txt.gz.tbi; done

#in case we later decide to VT decompose and normalize, annotate with VEP

#vt decompose $DATA/gnomad.exomes.r2.0.1.sites.vcf.gz -o /tmp/gnomad.exomes.r2.0.1.sites.dec.vcf
#vt normalize $DATA/gnomad.exomes.r2.0.1.sites.dec.vcf -o $DATA/gnomad.exomes.r2.0.1.sites.vt.vcf -r $DATA/grch37.fa
#perl $HOME/software/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl -i $DATA/gnomad.exomes.r2.0.1.sites.vt.vcf --cache --sift b --polyphen b --symbol --numbers --biotype --total_length --allele_number -o $DATA/gnomad.exomes.r2.0.1.sites.vep.vt.vcf --vcf --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,ALLELE_NUM --offline --fork 8
#get coverage and tbis from ExAC
