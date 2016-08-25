#cutting out introns for Monte Carlo analysis and uniqing in case of duplicates and sorting anyway:
bedtools intersect -a - -b <(zgrep protein_coding data/Homo_sapiens.GRCh37.75.gtf.gz | grep -Pw 'CDS|stop_codon' | bedtools merge) | sort -k1,1 -k2,2n | uniq | bedtools groupby -i stdin -g 1,2,3 -c 4,5,6,7,8,9,10,11,12,13 -o distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct | sort -k12,12nr
#previously matched on canonical transcript, in case this changes again, code is below
#bedtools intersect -a - -b $DATA/vepcanonicalexons.gtf -wb | awk 'match($29, /\"(\w*)/, t) {if ($4 == t[1]) print}' | cut -f -13 | sort -k1,1 -k2,2n | uniq
