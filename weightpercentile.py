import sys

totlen=0.0
f=open(sys.argv[1], "r")
f.readline()
for line in f:
    fields = line.strip().split("\t")
    totlen+=int(fields[2])-int(fields[1])
f.seek(0)
line = f.readline().strip() + "\t" + "weighted_pct"
pct=100.0; regionlength=0
print line 
for line in f:
    fields = line.strip().split("\t")
    regionlength += int(fields[2])-int(fields[1]) # a region can be split across multiple exons
    try:
        if fields[11]!=opct: # opct is original percentile, if it changes, we know the region is different and the region has changed and pct needs to be adjusted
            pct-=regionlength/totlen*100 # maybe the calculation should be from the bottom up?  previous calculation was.
            regionlength=0
    except NameError:
        opct = fields[11]
    opct = fields[11]
    line = line.strip() + "\t" + str(pct)
    print line
