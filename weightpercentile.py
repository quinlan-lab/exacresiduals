import sys
import toolshed as ts

totlen=0.0
for d in ts.reader(sys.argv[1]):
    totlen+=int(d['end'])-int(d['start']) 
f=open(sys.argv[1], "r")
line = ts.header(sys.argv[1])
line = "\t".join(line) + "\t" + "weighted_pct"
pct=100.0; regionlength=0
print line
for d in ts.reader(sys.argv[1], header='ordered'):
    regionlength += int(d['end'])-int(d['start'])
    try:
        if d['resid_pctile']!=opct:
            pct-=regionlength/totlen*100
            regionlength=0
    except NameError:
        opct=d['resid_pctile']
    opct=d['resid_pctile']
    line = "\t".join([i for i in d.values()]) + "\t" + str(pct)
    print line
