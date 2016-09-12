#creates flattened exome
grep -Pw 'stop_codon|CDS' /scratch/ucgd/lustre/u1021864/serial/Homo_sapiens.GRCh37.75.gtf | grep -w 'protein_coding' | sort -k1,1 -k4,4n | bedtools merge > flatexome.bed 

awk '{t+=$3-$2} END {print t}' flatexome.bed # length of flattened exome

awk '{t+=$3-$2} END {print t}' results/weightedresiduals.txt # length of final regions

#how much of the exome is covered by <10 bp regions
bedtools intersect -a <(awk '$3-$2<10' /scratch/ucgd/lustre/u1021864/serial/exac-regions.txt | cut -f -3) -b flatexome.bed | sort -k1,1 -k2,2n | bedtools merge | awk '{t+=$3-$2} END {print t}'

# total exome covered by our unfiltered regions (cut takes care of the 8 or 10 single bp regions missing fields for whatever reason)
bedtools intersect -a <(cut -f -3 /scratch/ucgd/lustre/u1021864/serial/exac-regions.txt) -b flatexome.bed  | sort -k1,1 -k2,2n | bedtools merge | awk '{t+=$3-$2} END {print t}'
