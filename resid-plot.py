import toolshed as ts
import matplotlib
matplotlib.use('Agg')
from matplotlib import gridspec
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
import cPickle as pickle
X = {"CpG": [], "gerp": []}

ys, genes = [], []
for i, d in enumerate(ts.reader(1)):
    try:
        gerps = [float(x) for x in d['gerp'].split(",")]
    except KeyError:
        continue
    genes.append((d['chrom'], str(d['start']), str(d['end']), d['gene'], d['transcript'], d['exon'], str(len(gerps))))
    
    coverage = map(float, d['coverage'].split(","))
    X['CpG'].append(float(d['cg_content']))
    X['gerp'].append(np.mean(gerps))

    ys.append(np.log(1.0 + np.sum(coverage)))
    #ys.append(np.sum(coverage))

gerp = X['gerp']
X['intercept'] = np.ones(len(ys))
del X['gerp']
X = pd.DataFrame(X)
#X = np.array(cgs)

results = sm.OLS(ys, X, hasconst=True).fit()
resid = OLSInfluence(results).get_resid_studentized_external()
variables={}
variables['cpg']=X['CpG']
variables['cov']=ys
variables['resid']=resid 
variables['rawresid']=results.resid
variables['genes']=genes
variables['gerp']=gerp
variables['intercept']=results.params['intercept']
variables['cpgcoef']=results.params['CpG']
pickle.dump(variables, open("var.pickle", "wb"))

resid_pctile = 100.0 * np.sort(resid).searchsorted(resid) / float(len(resid))
cov_pctile = 100.0 * np.sort(ys).searchsorted(ys) / float(len(ys))

print "chrom\tstart\tend\tgene\ttranscript\texon\tn\tcov_score\tcpg\tgerp_mean\tcov_cpg_resid\tcov_cpg_resid_pctile\tcov_pctile"
for i, row in enumerate(genes):
    print "\t".join(list(row) + ["%.3f" % ys[i], "%.3f" % X['CpG'][i], "%.3f" % gerp[i], "%.3f" % resid[i], "%.9f" % resid_pctile[i], "%.9f" % cov_pctile[i]
                                 ])
