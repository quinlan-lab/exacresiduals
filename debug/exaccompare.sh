#sed '1d' ../results/exacv1final/weightedresiduals-cpg-novariant.txt | sort -k1,1 -k2,2n | bgzip -c > exacccr.bed.gz; tabix exacccr.bed.gz
#sed '1d' ../results/parallel/weightedresiduals-cpg-novariant.txt | sort -k1,1 -k2,2n | bgzip -c > gnomadccr.bed.gz; tabix gnomadccr.bed.gz

bedtools intersect -a <(sed '1d' ../results/exacv1final/weightedresiduals-cpg-novariant.txt) -b ../../patho-exac.vcf.gz > exacpvars
bedtools intersect -a <(sed '1d' ../results/parallel/weightedresiduals-cpg-novariant.txt) -b ../../patho-gnomad.vcf.gz > gnomadpvars
bedtools intersect -a exacpvars -b gnomadpvars -wa -wb | cut -f 1,2,3,4,14,15,16,17,18,28 | awk '{print $0, $5-$10}' OFS='\t' | awk '$11 > 25' | sort -k11,11nr > pathodiff.txt
bedtools intersect -a <(sed '1d' ../results/exacv1final/weightedresiduals-cpg-novariant.txt) -b ../../benign-exac.vcf.gz > exacbvars
bedtools intersect -a <(sed '1d' ../results/parallel/weightedresiduals-cpg-novariant.txt) -b ../../benign-gnomad.vcf.gz > gnomadbvars
bedtools intersect -a exacbvars -b gnomadbvars -wa -wb | cut -f 1,2,3,4,14,15,16,17,18,28 | awk '{print $0, $10-$5}' OFS='\t' | awk '$11 > 25' | sort -k11,11nr > benigndiff.txt
python exaccompare.py
