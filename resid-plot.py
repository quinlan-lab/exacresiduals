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
from sklearn import preprocessing

import csv
import sys
csv.field_size_limit(34365000)
import cPickle as pickle
from cyvcf2 import VCF
import utils as u
from collections import defaultdict
X = defaultdict(list)
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cpg", help="cpg added to regression model", action="store_true", default=False)
parser.add_argument("-s", "--synonymous", help="synonymous added to regression model", action="store_true", default=False)
parser.add_argument("-f", "--file", help="regions input file, from exac-regions.py", required=True)
parser.add_argument("-n", "--nosingletons", help="if you do NOT want singletons", action="store_true", default=False)
parser.add_argument("-w", "--varflag", help="if you want separation by variant flags", action="store_true", default=False)

args=parser.parse_args()
cpg=args.cpg
synonymous=args.synonymous
nosingletons=args.nosingletons
rfile=args.file
varflag=args.varflag

exac=VCF('data/ExAC.r0.3.sites.vt.vep.vcf.gz') # only used for synonymous density calculation...can update with gnomAD if we really need it later
kcsq = exac["CSQ"]["Description"].split(":")[1].strip(' "').split("|")

ys, genes = [], []

def syn_density(pairs, d, exac, kcsq, nosingletons, varflag):
    syn=0
    synbool=False; prevvar=None
    if varflag:
        if 'VARTRUE' in d['varflag']: # don't need syn for a 0 bp region, i.e., variant, so give it the lowest possible, 0
            return syn
    for pair in pairs:
        if varflag:
            r0=str(int(pair[0])+1); r1=str(int(pair[1])); #in this case, does not include a variant at the end coordinate
        else:
            r0=str(int(pair[0])+1); r1=str(int(pair[1])-1);
        if not varflag:
            if int(r0)-int(r1)==1: continue # don't need syn for a region of length 1 (0 bp region), which it would be if a variant was included at the end coordinate
        for v in exac(d['chrom']+':'+r0+'-'+r1):
            if v.INFO['AC_Adj']==1 and nosingletons: continue
            if prevvar is not None and v.start+v.end==prevvar: continue
            if not (v.FILTER is None or v.FILTER == "PASS"): continue
            info = v.INFO
            try:
                csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]
            except KeyError:
                continue
            for csq in (c for c in csqs if c['BIOTYPE'] == 'protein_coding'):
                if csq['Feature'] == '' or csq['EXON'] == '' or csq['cDNA_position'] == '' or csq['SYMBOL']!=d['gene']: continue #in case non-exonic or not the same gene
                if not u.isfunctional(csq):
                    if not synbool:
                        syn+=1; synbool=True
                else:
                    if synbool:
                        syn-=1; break
            synbool=False
            prevvar=v.start+v.end

    return syn

varrow = []

for i, d in enumerate(ts.reader(rfile)):
    if d['chrom'] == 'X' or d['chrom'] == 'Y': continue
    pairs = [x.split("-") for x in d['ranges'].strip().split(",")]
    #try:
    #    if sum(e - s for s, e in (map(int, p) for p in pairs)) <= 10: 
    #        continue
    #except:
    #    print >>sys.stderr, d, pairs
    #    raise
    if 'VARTRUE' in d['varflag']:
        varrow.append((d['chrom'], str(d['start']), str(d['end']), d['gene'], d['transcript'], d['exon'], d['ranges'], d['varflag'], 0, 0))
        continue
    row=(d['chrom'], str(d['start']), str(d['end']), d['gene'], d['transcript'], d['exon'], d['ranges'], d['varflag'])
    if synonymous:
        syn=syn_density(pairs, d, exac, kcsq, nosingletons, varflag)
        if int(d['n_bases'])>1:
            if varflag:
                if 'VARTRUE' not in d['varflag']: # code here in case we decided to downweight differently later
                    d['syn_density']=syn/(float(d['n_bases'])-1); #+","+str(syn)+"/"+d['n_bases']; # -1 because we can't count the end coordinate, which is by default a variant
            else:
                d['syn_density']=syn/(float(d['n_bases'])-1); #+","+str(syn)+"/"+d['n_bases']; # -1 because we can't count the end coordinate, which is by default a variant
        else:
            if varflag:
                 if 'VARTRUE' not in d['varflag']: # code here in case we decided to downweight differently later
                    d['syn_density']=0
            else:
                d['syn_density']=0
        X['syn'].append(float(d['syn_density'])) # 1-syn if we want to use as a measure of constraint; syn as a measure of mutability
        row = row + ("%.3f" % float(d['syn_density']),)
    else:
        d['syn_density']="na" # calculating synonymous density is really slow, so if we don't need to, we'd rather not.
        row = row + (d['syn_density'],)
    if cpg:
        if varflag:
            if 'VARTRUE' not in d['varflag']: # code here in case we decided to downweight differently later
                X['CpG'].append(float(d['cg_content']))
        else: 
            X['CpG'].append(float(d['cg_content']))
    row = row + ("%.3f" % float(d['cg_content']),)
    genes.append(row)
    coverage=[]
    for val in d['coverage'].split(","):
        if val:
            if varflag:
                if 'VARTRUE' not in d['varflag']: # code here in case we decided to downweight differently later
                    coverage.append(float(val))
            else:
                coverage.append(float(val))
    if not coverage:
        if varflag:
             if 'VARTRUE' not in d['varflag']: # code here in case we decided to downweight differently later
                coverage=[0]
        else:
            coverage=[0]
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

lowestresidual=np.min(resid)-.001
#for i, row in enumerate(genes):
#    if "VARTRUE" in row[7] and varflag: #row[7] is varflag
#        resid[i]=lowestresidual
resid=resid.tolist()
for i, row in enumerate(varrow):
    resid.append(lowestresidual)
    genes.append(row)
    ys.append(0)
    
X_train=np.array(resid).reshape(len(resid),1)
min_max_scaler = preprocessing.MinMaxScaler(feature_range=(0,100))
resid_pctile = min_max_scaler.fit_transform(X_train)
#resid_pctile = 100.0 * np.sort(resid).searchsorted(resid) / float(len(resid))

assert len(genes) == len(ys) == len(resid)

print "chrom\tstart\tend\tgene\ttranscript\texon\tranges\tvarflag\tsyn_density\tcpg\tcov_score\tresid\tresid_pctile"
for i, row in enumerate(genes):
    #if "VARTRUE" in row[7] and varflag: #row[7] is varflag 
    vals = ["%.3f" % ys[i], "%.3f" % resid[i], "%.9f" % resid_pctile[i]]
    #if not "," in row[-1]:
    #    if not row[-1]:
    #        row=list(row)
    #        row[-1]=row[1]+"-"+row[2]
    #    print "\t".join(list(row) + vals)
    #    continue
     
    ranges = [x.split("-") for x in row[6].split(",")]
    row=list(row)
    for s, e in ranges:
        row[1], row[2] = s, e
        print "\t".join(map(str,list(row) + vals))
