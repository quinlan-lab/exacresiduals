# SEE: https://github.com/quinlan-lab/lab-wiki/blob/master/projects/residuals.md
import sys

from collections import defaultdict
from cyvcf2 import VCF
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import OLSInfluence

exac = VCF('/usr/local/src/gemini_install/data/gemini_data/ExAC.r0.3.sites.vep.tidy.vcf.gz')

# CSQ keys
kcsq = exac["CSQ"]["Description"].split(":")[1].strip(' "').split("|")

exons = defaultdict(list)


def isfunctional(csq):
    return any(c in ('stop_gained', 'stop_lost', 'start_lost', 'initiator_codon_variant', 'rare_amino_acid_variant',
                     'missense_variant', '5_prime_UTR_premature_start_codon_gain_variant', 'protein_altering_variant',
                     'splice_acceptor_variant', 'splice_donor_variant', 'frameshift_variant')
               for c in csq['Consequence'].split('&'))

def exonkey(csq):
    if csq['Feature'] == '':
        return None
    if csq['EXON'] == '':
        return None
        if csq['INTRON'] == '':
            return None
        else:
            return "intron"+csq['INTRON'], csq['Feature']

    return csq['Feature'], csq['EXON']

raise Exception("TODO: get GTF of exons and add that to output along with re-positioned variant within the transcripts")
raise Exception("used bed files of coverage (read into InterLap) to limit to 80% of samples > 10X")
print "#chrom\tstart\tend\taf\tfunctional\tgene\ttranscript\texon\timpact"
for i, v in enumerate(exac, start=1):
    info = v.INFO
    try:
        csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]
    except KeyError:
        continue
    af = info['AC_Adj'] / float(info['AN_Adj'] or 1)
    try:
        csq = next(c for c in csqs if c['CANONICAL'] == 'YES' and c['Allele'] == v.ALT[0])
    except StopIteration:
        continue

    key = exonkey(csq)
    if key is None: continue

    print "\t".join(map(str, (v.CHROM, v.start, v.end, af, int(isfunctional(csq)),
        csq['SYMBOL'], csq['Feature'], csq['EXON'], csq['Consequence'])))
