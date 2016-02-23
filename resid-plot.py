import toolshed as ts
from matplotlib import pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import OLSInfluence

import csv
import sys
csv.field_size_limit(sys.maxsize)


def fn(scores, covs):
    return sum(sc * cov for sc, cov in zip(scores, covs))

xs, ys, genes = [], [], []
for i, d in enumerate(ts.reader(1)):
    genes.append((d['chrom'], str(d['start']), str(d['end']), d['gene'], d['exon']))

    gerps = map(float, d['gerp'].split(","))
    coverage = map(float, d['coverage'].split(","))

    score = fn(gerps, coverage)

    ys.append(score)
    xs.append(sum(coverage))
    if i >= 3089770: break

results = sm.OLS(ys, xs, hasconst=False).fit()

l = results.predict([0, max(xs)])

plt.plot(xs, ys, marker='.', ls='none')
plt.plot([0, max(xs)], l)
plt.title("split by missense")
plt.xlabel("sum(coverage)")
plt.ylabel("sum(cov * g for cov, g in zip(coverage, gerp))")
plt.show()

resid = results.resid

for i, row in enumerate(genes):
    print "\t".join(list(row) + ["%.3f" % resid[i]])
