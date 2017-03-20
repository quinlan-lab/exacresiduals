#!/bin/bash
#SBATCH --account=quinlan-kp
#SBATCH --partition=quinlan-kp
#SBATCH -o ../%j-%N.out
#SBATCH -e ../%j-%N.err
#SBATCH --time=6:00:00
set -exo pipefail -o nounset

## folder made by date, in case we make major changes to exac-regions.py or resid-plot.py ##
date=2017_03_16
#mkdir -p results/$date/
## generates regions and residuals files ##
python exac-regions.py > results/$date/exac-regions.txt
python resid-plot.py results/$date/exac-regions.txt > results/$date/resids.txt
cat <(head -1 results/$date/resids.txt) <(sed '1d' results/$date/resids.txt | sort -k12,12nr) > /tmp/residsort.txt
python weightpercentile.py /tmp/residsort.txt > results/$date/weightedresiduals.txt
#no-singletons
#python exac-regions.py -s > results/$date/exac-regions-nosingletons.txt
#python resid-plot.py results/$date/exac-regions-nosingletons.txt > results/$date/resids-nosingletons.txt
#cat <(head -1 results/$date/resids-nosingletons.txt) <(sed '1d' results/$date/resids-nosingletons.txt | sort -k12,12nr) > /tmp/residsort-nosingletons.txt
#python weightpercentile.py /tmp/residsort-nosingletons.txt > results/$date/weightedresiduals-nosingletons.txt
#getting unfiltered regions, purely exonic (comment out coverage, self-chain, and seg dup filters in utils.py)
#python exac-regions.py > results/$date/unfilteredregions.txt
#python resid-plot.py results/$date/unfilteredregions.txt > results/$date/unfilteredresiduals.txt
## getting exonic-only residuals and getting top residuals by percentile and middle residual regions by exonic BP totals and closeness to 0 raw resid values ##
#sed '1d' results/weightedresiduals.txt | awk '$14 >= 99' > ../regions/topresid.txt
#python middle.py -b > ../regions/midresid.txt # -b to run by total basepair matching at the ~0 residual score line (default); -g to match by number of genes for gene comparison

## old code for bottom residuals and genewide stuff, may need some editing ##
#sed '1d' results/$date/exonicresiduals.txt | awk '$12 <=1' > ../regions/topresid.txt
#sed '1d' results/$date/exonicresiduals.txt | sort -k4,4 | bedtools groupby -g 4 -c 1,2,3,8 -o distinct,min,max,mean | awk '{print $2,$3,$4,$1,$5}' OFS="\t" | awk '{gene[$4]=$5; row[$4]=$0} END {for (i in gene) for (j in gene) {if (gene[i]>=gene[j]) rank[i]+=1}} END {for (i in rank) print row[i],100*rank[i]/NR}' > meangeneresiduals.txt
#sed '1d' results/$date/exonicresiduals.txt | sort -k4,4 | bedtools groupby -g 4 -c 1,2,3,8 -o distinct,min,max,max | awk '{print $2,$3,$4,$1,$5}' OFS="\t" | awk '{gene[$4]=$5; row[$4]=$0} END {for (i in gene) for (j in gene) {if (gene[i]>=gene[j]) rank[i]+=1}} END {for (i in rank) print row[i],100*rank[i]/NR}' > maxgeneresiduals.txt
#awk '$6 > 99' meangeneresiduals.txt | bedtools intersect -a stdin -b $DATA/vepcanonicalexons.gtf > ../regions/topmeangeneresid.txt
#awk '$5 > -0.033 && $5 < 0.033' meangeneresiduals.txt | bedtools intersect -a stdin -b $DATA/vepcanonicalexons.gtf > ../regions/midmeangeneresid.txt
#awk '$6 > 99' maxgeneresiduals.txt | bedtools intersect -a stdin -b $DATA/vepcanonicalexons.gtf > ../regions/topmaxgeneresid.txt
#awk '$6 < 5' maxgeneresiduals.txt | bedtools intersect -a stdin -b $DATA/vepcanonicalexons.gtf > ../regions/midmaxgeneresid.txt
