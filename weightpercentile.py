import sys
import toolshed as ts
from sklearn import preprocessing

regions, weighted = [], []
totlen=0.0

for d in ts.reader(sys.argv[1]):
    totlen+=int(d['end'])-int(d['start'])
    regions.append(d)

header = ts.header(sys.argv[1])
print "\t".join(header) + "\t" + "weighted_pct"
pct=100.0; regionlength=0

opct = regions[0]['resid_pctile']
for d in regions:
    regionlength += int(d['end'])-int(d['start'])
    if d['resid_pctile']!=opct:
        pct-=regionlength/totlen*100
        regionlength=0
        opct=d['resid_pctile']
    weighted.append(pct)
#    line = "\t".join([i for i in d.values()]) + "\t" + str(pct)
#    print line
min_max_scaler = preprocessing.MinMaxScaler()
resid_pctile = list(100*min_max_scaler.fit_transform(weighted))
for i, d in enumerate(regions):
    print "\t".join(d[h] for h in header) + "\t" + str(resid_pctile[i])
