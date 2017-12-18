sed '1d' results/gnomAD5x.1/weightedresiduals-cpg-novariant.txt | sort -k1,1 -k2,2n | bgzip -c -@ 12 > gnomad5x.1-ccrs.bed.gz; tabix -f gnomad5x.1-ccrs.bed.gz
sed '1d' results/gnomAD5x.9/weightedresiduals-cpg-novariant.txt | sort -k1,1 -k2,2n | bgzip -c -@ 12 > gnomad5x.9-ccrs.bed.gz; tabix -f gnomad5x.9-ccrs.bed.gz
sed '1d' results/gnomAD50x.1/weightedresiduals-cpg-novariant.txt | sort -k1,1 -k2,2n | bgzip -c -@ 12 > gnomad50x.1-ccrs.bed.gz; tabix -f gnomad50x.1-ccrs.bed.gz
sed '1d' results/gnomAD50x.9/weightedresiduals-cpg-novariant.txt | sort -k1,1 -k2,2n | bgzip -c -@ 12 > gnomad50x.9-ccrs.bed.gz; tabix -f gnomad50x.9-ccrs.bed.gz
sed '1d' results/gnomAD10x.5/weightedresiduals-cpg-novariant.txt | sort -k1,1 -k2,2n | bgzip -c -@ 12 > gnomad10x.5-ccrs.bed.gz; tabix -f gnomad10x.5-ccrs.bed.gz
sed '1d' results/gnomAD30x.5/weightedresiduals-cpg-novariant.txt | sort -k1,1 -k2,2n | bgzip -c -@ 12 > gnomad30x.5-ccrs.bed.gz; tabix -f gnomad30x.5-ccrs.bed.gz
