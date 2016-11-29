import toolshed as ts
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import OLSInfluence
import scipy.stats as ss
from statsmodels.formula.api import ols
import pandas as pd
from scipy.stats.mstats import hmean


import csv
import sys
csv.field_size_limit(14365000)
import cPickle as pickle
X = {"CpG": []}

ys, genes = [], []
for i, d in enumerate(ts.reader(1)):
    if d['chrom'] == 'X' or d['chrom'] == 'Y': continue
    if int(d['end']) - int(d['start']) < 10: continue

    pairs = [x.split("-") for x in d['ranges'].strip().split(",")]
    try:
        if sum(e - s for s, e in (map(int, p) for p in pairs)) <= 10: # = is because the ranges are in VCF space, not BED space
            continue
    except:
        print >>sys.stderr, d, pairs
        raise

    genes.append((d['chrom'], str(d['start']), str(d['end']), d['gene'], d['transcript'], d['exon'], d['ranges']))

    coverage = map(float, d['coverage'].split(","))
    X['CpG'].append(float(d['cg_content']))
    ys.append(sum(coverage))

X['intercept'] = np.ones(len(ys))
X = pd.DataFrame(X)


results = sm.OLS(ys, X, hasconst=True).fit()
resid = OLSInfluence(results).get_resid_studentized_external()
#variables={}
#variables['cpg']=X['CpG']
#variables['cov']=ys
#variables['resid']=resid
#variables['rawresid']=results.resid
#variables['genes']=genes
#variables['gerp']=gerp
#variables['intercept']=results.params['intercept']
#variables['cpgcoef']=results.params['CpG']
#pickle.dump(variables, open("var.pickle", "wb"))

resid_pctile = 100.0 * np.sort(resid).searchsorted(resid) / float(len(resid))

assert len(genes) == len(ys) == len(resid)

print "chrom\tstart\tend\tgene\ttranscript\texon\tranges\tcov_score\tcpg\tcov_cpg_resid\tcov_cpg_resid_pctile"
for i, row in enumerate(genes):
    vals = ["%.3f" % ys[i], "%.3f" % X['CpG'][i], "%.3f" % resid[i], "%.9f" % resid_pctile[i]]
    if not "," in row[-1]:
        print "\t".join(list(row) + vals)
        continue
     
    ranges = [x.split("-") for x in row[-1].split(",")]
    row=list(row)
    for s, e in ranges:
        row[1], row[2] = s, e
        print "\t".join(list(row) + vals)
