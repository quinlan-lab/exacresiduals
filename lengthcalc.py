import sys

for line in sys.stdin:
    fields = line.strip().split("\t")
    gene = fields[0]; ranges = fields[1]; resid = fields[2]
    le = [l.split("-") for l in ranges.split(",")]
    length = 0
    for range in le:
        length += int(range[-1]) - int(range[0])
    print "\t".join([gene, ranges, resid, str(length)]) 
