#!/bin/bash
#SBATCH --account=quinlan-kp
#SBATCH --partition=quinlan-kp
#SBATCH -o ../%j-%N.out
#SBATCH -e ../%j-%N.err
#SBATCH --time=24:00:00
#set -exo pipefail -o nounset

# get ExAC v2.0 VCF, coverage files
#wget -P $DATA https://storage.googleapis.com/gnomad-public/release-170228/vcf/exomes/gnomad.exomes.r2.0.1.sites.vcf.gz
#for i in {1..22} X Y; do wget -P $HOME/analysis/exacresiduals/data https://storage.googleapis.com/gnomad-public/release-170228/coverage/exomes/exacv2.chr$i.cov.txt.gz; done
#for i in {1..22} X Y; do wget -P $HOME/analysis/exacresiduals/data https://storage.googleapis.com/gnomad-public/release-170228/coverage/exomes/exacv2.chr$i.cov.txt.gz.tbi; done

#gnomAD

bash $HOME/analysis/varmake.sh $DATA/gnomad.exomes.r2.0.1.sites.vcf.gz

#vt decompose $DATA/gnomad.exomes.r2.0.1.sites.vcf.gz -o $DATA/gnomad.exomes.r2.0.1.sites.dec.vcf -s 
#vt normalize $DATA/gnomad.exomes.r2.0.1.sites.dec.vcf -o $DATA/gnomad.exomes.r2.0.1.sites.vt.vcf -r $DATA/grch37.fa
#perl $HOME/software/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl -i $DATA/gnomad.exomes.r2.0.1.sites.vt.vcf --cache --sift b --polyphen b --symbol --numbers --biotype --allele_number -o $DATA/gnomad.exomes.r2.0.1.sites.vep.vt.vcf --vcf --fields ALLELE,Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,ALLELE_NUM,cDNA_position --offline --fork 12 --force_overwrite
#bgzip $DATA/gnomad.exomes.r2.0.1.sites.vep.vt.vcf
#tabix $DATA/gnomad.exomes.r2.0.1.sites.vep.vt.vcf.gz

# old ExAC
vt decompose $DATA/ExAC.r1.sites.vep.vcf.gz -o $DATA/ExAC.r1.sites.dec.vcf -s
vt normalize $DATA/ExAC.r1.sites.dec.vcf -o $DATA/ExAC.r1.sites.vt.vcf -r $DATA/grch37.fa
perl $HOME/software/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl -i $DATA/ExAC.r1.sites.vt.vcf --cache --sift b --polyphen b --symbol --numbers --biotype --allele_number -o $DATA/ExAC.r1.sites.vt.vep.vcf --vcf --fields ALLELE,Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,ALLELE_NUM,cDNA_position --offline --fork 12 --force_overwrite
bgzip $DATA/ExAC.r1.sites.vt.vep.vcf
tabix $DATA/ExAC.r1.sites.vt.vep.vcf.gz
