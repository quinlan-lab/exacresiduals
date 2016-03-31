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


import csv
import sys
csv.field_size_limit(14365000)

X = {"CpG": [], "gerp": []}

ys, genes = [], []
for i, d in enumerate(ts.reader(1)):
    gerps = [float(x) for x in d['gerp'].split(",")]
    genes.append((d['chrom'], str(d['start']), str(d['end']), d['gene'], d['exon'], str(len(gerps))))


    coverage = map(float, d['coverage'].split(","))
    X['CpG'].append(float(d['cg_content']))
    X['gerp'].append(np.mean(gerps))

    ys.append(np.log(1.0 + np.sum(coverage)))

del X['gerp']
X = pd.DataFrame(X)
#X = np.array(cgs)

results = sm.OLS(ys, X, hasconst=False).fit()
resid = OLSInfluence(results).get_resid_studentized_external()

resid_pctile = 100.0 * np.sort(resid).searchsorted(resid) / float(len(resid))
cov_pctile = 100.0 * np.sort(ys).searchsorted(ys) / float(len(ys))

print "chrom\tstart\tend\tgene\texon\tn\tgerp_cpg_resid\tgerp_cpg_resid_pctile\tcov_pctile"
for i, row in enumerate(genes):
    print "\t".join(list(row) + ["%.3f" % resid[i], "%.9f" % resid_pctile[i], "%.9f" % cov_pctile[i]
                                 ])
