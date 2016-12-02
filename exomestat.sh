DATA=/scratch/ucgd/lustre/u1021864/serial/
#creates flattened transcriptome
grep -P 'protein_coding\tstop_codon|protein_coding\tCDS' $DATA/Homo_sapiens.GRCh37.75.gtf | grep -P "^1|^2|^3|^4|^5|^6|^7|^8|^9|^10|^11|^12|^13|^14|^15|^16|^17|^18|^19|^20|^21|^22" | awk '{$4=$4-1; print $0}' OFS='\t' | cut -f 1,4,5,12,16 | sort -k5,5 -k1,1 -k2,2n > codingtranscriptome.bed
python flattenexome.py | sort -k1,1 -k2,2n > flatexome.bed

echo "length of our flattened exome"
awk '{t+=$3-$2} END {print t}' flatexome.bed # length of flattened exome
echo "length of our most recent regions"
awk '{t+=$3-$2} END {print t}' results/2016_12_01/weightedresiduals.txt # length of final regions

cat coveragecut.txt segdupsremoved.txt selfchainremoved.txt | cut -f -3 | sort -k1,1 -k2,2n | bedtools merge | bedtools intersect -a - -b flatexome.bed -sorted > cutregions.txt
awk '{t+=$3-$2} END {print t}' cutregions.txt

bedtools intersect -a selfchainremoved.txt -b cutregions.txt | sort -k1,1 -k2,2n | awk '{t+=$3-$2} END {print t}'
bedtools intersect -a segdupsremoved.txt -b cutregions.txt | sort -k1,1 -k2,2n | awk '{t+=$3-$2} END {print t}'
bedtools intersect -a coveragecut.txt -b cutregions.txt | sort -k1,1 -k2,2n | awk '{t+=$3-$2} END {print t}'
# amount that overlaps with the other two:
bedtools intersect -a <(bedtools intersect -a coveragecut.txt -b flatexome.bed) -b selfchainremoved.txt segdupsremoved.txt | sort -k1,1 -k2,2n | bedtools merge | awk '{t+=$3-$2} END {print t}'
bedtools intersect -a <(bedtools intersect -a segdupsremoved.txt -b flatexome.bed) -b coveragecut.txt selfchainremoved.txt | sort -k1,1 -k2,2n | bedtools merge | awk '{t+=$3-$2} END {print t}'
bedtools intersect -a <(bedtools intersect -a selfchainremoved.txt -b flatexome.bed) -b coveragecut.txt segdupsremoved.txt | sort -k1,1 -k2,2n | bedtools merge | awk '{t+=$3-$2} END {print t}'
