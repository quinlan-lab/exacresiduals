#!/bin/bash
#SBATCH --account=quinlan-kp
#SBATCH --partition=quinlan-kp
#SBATCH -o ../%j-%N.out
#SBATCH -e ../%j-%N.err
#SBATCH --time=4:00:00
set -exo pipefail -o nounset

#folder made by date, in case we make major changes to exac-regions.py or resid-plot.py
date=2016_06_15
mkdir -p results/$date/
#generates regions and residuals files
python exac-regions.py | bgzip -c > results/$date/regs.bed.gz
python resid-plot.py results/$date/regs.bed.gz > results/$date/resids.txt
#getting exonic-only residuals and getting top residuals by percentile and middle residual regions by exonic BP totals and closeness to 0 raw resid values
cat <(head -1 results/$date/resids.txt) <(sed '1d' results/$date/resids.txt | bash dups.sh) > results/$date/exonicresiduals.txt
sed '1d' results/$date/exonicresiduals.txt | awk '$12 >= 99' > ../regions/topresid.txt
python middle.py > ../regions/midresid.txt

# old code for bottom residuals and genewide stuff, may need some editing
#sed '1d' results/$date/exonicresiduals.txt | awk '$12 <=1' > ../regions/topresid.txt
#sed '1d' results/$date/exonicresiduals.txt | sort -k4,4 | bedtools groupby -g 4 -c 1,2,3,8 -o distinct,min,max,mean | awk '{print $2,$3,$4,$1,$5}' OFS="\t" | awk '{gene[$4]=$5; row[$4]=$0} END {for (i in gene) for (j in gene) {if (gene[i]>=gene[j]) rank[i]+=1}} END {for (i in rank) print row[i],100*rank[i]/NR}' > meangeneresiduals.txt
#sed '1d' results/$date/exonicresiduals.txt | sort -k4,4 | bedtools groupby -g 4 -c 1,2,3,8 -o distinct,min,max,max | awk '{print $2,$3,$4,$1,$5}' OFS="\t" | awk '{gene[$4]=$5; row[$4]=$0} END {for (i in gene) for (j in gene) {if (gene[i]>=gene[j]) rank[i]+=1}} END {for (i in rank) print row[i],100*rank[i]/NR}' > maxgeneresiduals.txt
#awk '$6 > 99' meangeneresiduals.txt | bedtools intersect -a stdin -b $DATA/vepcanonicalexons.gtf > ../regions/topmeangeneresid.txt
#awk '$5 > -0.033 && $5 < 0.033' meangeneresiduals.txt | bedtools intersect -a stdin -b $DATA/vepcanonicalexons.gtf > ../regions/midmeangeneresid.txt
#awk '$6 > 99' maxgeneresiduals.txt | bedtools intersect -a stdin -b $DATA/vepcanonicalexons.gtf > ../regions/topmaxgeneresid.txt
#awk '$6 < 5' maxgeneresiduals.txt | bedtools intersect -a stdin -b $DATA/vepcanonicalexons.gtf > ../regions/midmaxgeneresid.txt
