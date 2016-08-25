import sys

totlen=0.0
f=open(sys.argv[1], "r")
f.readline()
for line in f:
    fields = line.strip().split("\t")
    totlen+=int(fields[2])-int(fields[1])
f.seek(0)
line = f.readline().strip() + "\t" + "weighted_pct"
pct=100.0
print line 
for line in f:
    fields = line.strip().split("\t")
    regionlength = int(fields[2])-int(fields[1])
    try:
        if fields[11]!=opct:
            pct=pct-regionlength/totlen
    except NameError:
        opct = fields[11]
    opct = fields[11]
    line = line.strip() + "\t" + str(pct)
    print line
