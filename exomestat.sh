grep -Pw 'stop_codon|CDS' /scratch/ucgd/lustre/u1021864/serial/Homo_sapiens.GRCh37.75.gtf | grep -w 'protein_coding' | sort -k1,1 -k4,4n | bedtools merge | awk '{t+=$3-$2} END {print t}'

awk '{t+=$3-$2} END {print t}' results/weightedresiduals.txt
