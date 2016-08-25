#To get residual files here use branch test:

Use https://github.com/quinlan-lab/exacresiduals/tree/test (branch test is the most up-to-date)

Folder located at: /uufs/chpc.utah.edu/common/home/u1021864/analysis/exacresiduals ("exacresiduals")

##Under the folder "exacresiduals":

The folder "../regions" (this folder should be in the folder "analysis" from https://github.com/quinlan-lab/regionanalysis) is where the script "regions.sh" creates and stores the "topresid.txt" and "midresid.txt" files for Monte Carlo analysis. At one point we were trying to create a gene based mean score for residuals and those are stored in there for now as well.  The script now accesses "weightpercentile.py" which is used to create the size-weighted percentiles for the new residual file which is stored in the folder "results" as "weightedresiduals.txt".

The script "dups.sh" flattens the regions so each bp in the exome is only represented once and uses bedtools groupby to ensure that genes with different percentile scores for the same exome-space region are represented only once as well.  This script is accessed by "regions.sh" to create both the "exonicresiduals.txt" file in the "results/$date" folder put also the "weightedresiduals.txt" by consequence of its use by "weightpercentile.py".

