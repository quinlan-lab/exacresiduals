#cutting out introns for Monte Carlo analysis and uniqing in case of duplicates and sorting anyway:
bedtools intersect -a - -b <(zgrep protein_coding data/Homo_sapiens.GRCh37.75.gtf.gz | grep -Pw 'CDS|stop_codon') | sort -k1,1 -k2,2n | uniq
#previously matched on canonical transcript, in case this changes again, code is below
#bedtools intersect -a - -b $DATA/vepcanonicalexons.gtf -wb | awk 'match($29, /\"(\w*)/, t) {if ($4 == t[1]) print}' | cut -f -13 | sort -k1,1 -k2,2n | uniq
