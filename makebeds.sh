sed '1d' results/gnomAD1x.1syn/weightedresiduals-cpg-synonymous-novariant.txt | sort -k1,1 -k2,2n | bgzip -c -@ 12 > gnomad1x.1syn-ccrs.bed.gz; tabix -f gnomad1x.1syn-ccrs.bed.gz
sed '1d' results/gnomAD1x.9syn/weightedresiduals-cpg-synonymous-novariant.txt | sort -k1,1 -k2,2n | bgzip -c -@ 12 > gnomad1x.9syn-ccrs.bed.gz; tabix -f gnomad1x.9syn-ccrs.bed.gz
sed '1d' results/gnomAD5x.1syn/weightedresiduals-cpg-synonymous-novariant.txt | sort -k1,1 -k2,2n | bgzip -c -@ 12 > gnomad5x.1syn-ccrs.bed.gz; tabix -f gnomad5x.1syn-ccrs.bed.gz
sed '1d' results/gnomAD5x.9syn/weightedresiduals-cpg-synonymous-novariant.txt | sort -k1,1 -k2,2n | bgzip -c -@ 12 > gnomad5x.9syn-ccrs.bed.gz; tabix -f gnomad5x.9syn-ccrs.bed.gz
sed '1d' results/gnomAD50x.1syn/weightedresiduals-cpg-synonymous-novariant.txt | sort -k1,1 -k2,2n | bgzip -c -@ 12 > gnomad50x.1syn-ccrs.bed.gz; tabix -f gnomad50x.1syn-ccrs.bed.gz
sed '1d' results/gnomAD50x.9syn/weightedresiduals-cpg-synonymous-novariant.txt | sort -k1,1 -k2,2n | bgzip -c -@ 12 > gnomad50x.9syn-ccrs.bed.gz; tabix -f gnomad50x.9syn-ccrs.bed.gz
sed '1d' results/gnomAD30x.5syn/weightedresiduals-cpg-synonymous-novariant.txt | sort -k1,1 -k2,2n | bgzip -c -@ 12 > gnomad30x.5syn-ccrs.bed.gz; tabix -f gnomad30x.5syn-ccrs.bed.gz
sed '1d' results/gnomAD10x.5syn/weightedresiduals-cpg-synonymous-novariant.txt | sort -k1,1 -k2,2n | bgzip -c -@ 12 > gnomad10x.5syn-ccrs.bed.gz; tabix -f gnomad10x.5syn-ccrs.bed.gz
sed '1d' results/Xchromonly/weightedresiduals-cpg-synonymous-novariant.txt | sort -k1,1 -k2,2n | bgzip -c -@ 12 > xchrom-ccrs.bed.gz; tabix -f xchrom-ccrs.bed.gz
#awk 'NR==1{printf "#"} {print}' results/gnomAD10x.5syn/weightedresiduals-cpg-synonymous-novariant.txt | sort -k1,1 -k2,2n | bgzip -c -@ 12 > gnomad10x.5syn-ccrs.bed.gz; tabix -f gnomad10x.5syn-ccrs.bed.gz
sed '1d' results/ExACv1syn/weightedresiduals-cpg-synonymous-novariant.txt | sort -k1,1 -k2,2n | bgzip -c -@ 12 > exacv1syn-ccrs.bed.gz; tabix -f exacv1syn-ccrs.bed.gz
#awk 'NR==1{printf "#"} {print}' results/ExACv1syn/weightedresiduals-cpg-synonymous-novariant.txt | sort -k1,1 -k2,2n | bgzip -c -@ 12 > exacv1syn-ccrs.bed.gz; tabix -f exacv1syn-ccrs.bed.gz
