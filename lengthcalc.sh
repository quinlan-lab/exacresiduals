zcat gnomad10x.5syn-ccrs.bed.gz | cut -f 4,7,14 | python lengthcalc.py > lengths
