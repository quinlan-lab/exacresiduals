sort -k11,11nr results/2016_06_16/resids.txt | head -100 | cut -f -3 | awk '{print $1":"$2"-"$3}' > hundred.txt    
for i in $(cat hundred.txt); do    
    python plot-cov.py $i
done
python makepdf.py hundred.txt
