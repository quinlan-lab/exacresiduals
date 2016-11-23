DATA=/scratch/ucgd/lustre/u1021864/serial/
#creates flattened exome
grep -Pw 'stop_codon|CDS' $DATA/Homo_sapiens.GRCh37.75.gtf | grep -w 'protein_coding' | sort -k1,1 -k4,4n | bedtools merge | awk '{print $1"\t"$2-1"\t"$3}'| grep -P "^1|^2|^3|^4|^5|^6|^7|^8|^9|^10|^11|^12|^13|^14|^15|^16|^17|^18|^19|^20|^21|^22" > flatexome.bed # removed X and Y chromosomes since we can't model them at the moment also removed MT, HGPATCH, etc.

echo "length of our flattened exome"
awk '{t+=$3-$2} END {print t}' flatexome.bed # length of flattened exome
echo "length of our most recent regions"
awk '{t+=$3-$2} END {print t}' results/2016_11_17/weightedresiduals.txt # length of final regions

# total exome covered by our unfiltered regions (cut takes care of the 8 or 10 single bp regions missing fields for whatever reason)
bedtools intersect -a <(cut -f -3 $DATA/exac-regions.txt | grep -Pv '^X|^Y') -b flatexome.bed  | sort -k1,1 -k2,2n | bedtools merge | awk '{t+=$3-$2} END {print t}'

cat coveragecut.txt segdupsremoved.txt selfchainremoved.txt | cut -f -3 | bedtools intersect -a - -b flatexome.bed | sort -k1,1 -k2,2n | bedtools merge > cutregions.txt
awk '{t+=$3-$2} END {print t}' cutregions.txt

bedtools intersect -a selfchainremoved.txt -b cutregions.txt | sort -k1,1 -k2,2n | bedtools merge | awk '{t+=$3-$2} END {print t}'
bedtools intersect -a segdupsremoved.txt -b cutregions.txt | sort -k1,1 -k2,2n | bedtools merge | awk '{t+=$3-$2} END {print t}'
bedtools intersect -a coveragecut.txt -b cutregions.txt | sort -k1,1 -k2,2n | bedtools merge | awk '{t+=$3-$2} END {print t}'
# amount that overlaps with the other two:
bedtools intersect -a <(bedtools intersect -a coveragecut.txt -b flatexome.bed) -b selfchainremoved.txt segdupsremoved.txt | sort -k1,1 -k2,2n | bedtools merge | awk '{t+=$3-$2} END {print t}'
bedtools intersect -a <(bedtools intersect -a segdupsremoved.txt -b flatexome.bed) -b coveragecut.txt selfchainremoved.txt | sort -k1,1 -k2,2n | bedtools merge | awk '{t+=$3-$2} END {print t}'
bedtools intersect -a <(bedtools intersect -a selfchainremoved.txt -b flatexome.bed) -b coveragecut.txt segdupsremoved.txt | sort -k1,1 -k2,2n | bedtools merge | awk '{t+=$3-$2} END {print t}'
#amount that is removed from filtering <10 bp
awk '{t+=$3-$2} END {print t}' results/2016_11_17/weightedresiduals.txt
awk '{t+=$3-$2} END {print t}' results/2016_11_17/nosizefilter/weightedresiduals.txt
