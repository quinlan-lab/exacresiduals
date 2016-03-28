import toolshed as ts
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

cutoff = 0.3

ys, genes = [], []
X = {'CpG': [], "gerp": []}
for i, d in enumerate(ts.reader(1)):
        gerps = [float(x) for x in d['gerp'].split(",")]

        genes.append((d['chrom'], str(d['start']), str(d['end']), d['gene'],
            d['exon'], str(len(gerps))))

        coverage = map(float, d['coverage'].split(","))
        ys.append(sum(c for c in coverage if c > cutoff))
        X['CpG'].append(float(d['cg_content']))
        X['gerp'].append(np.mean(gerps))


X = pd.DataFrame(X)
results = sm.OLS(ys, X, hasconst=False).fit()

print >>sys.stderr, results.params

resid = results.resid
resid = OLSInfluence(results).get_resid_studentized_external()

pctile = 100.0 * np.sort(resid).searchsorted(resid) / float(len(resid))

cov_score_pctile = 100.0 * np.sort(ys).searchsorted(ys) / float(len(ys))

print "chrom\tstart\tend\tgene\texon\tn\tgerp_cpg_resid\tgerp_cpg_resid_pctile\tcov_score_pctile"
for i, row in enumerate(genes):
    print "\t".join(list(row) + ["%.3f" % resid[i], "%.9f" % pctile[i],
                                 "%.9f" % cov_score_pctile[i],
                                 ])

