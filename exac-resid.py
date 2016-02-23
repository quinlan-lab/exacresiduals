# SEE: https://github.com/quinlan-lab/lab-wiki/blob/master/projects/residuals.md
import sys
from collections import defaultdict

from cyvcf2 import VCF

exac = VCF('/usr/local/src/gemini_install/data/gemini_data/ExAC.r0.3.sites.vep.tidy.vcf.gz')

# CSQ keys
kcsq = exac["CSQ"]["Description"].split(":")[1].strip(' "').split("|")

exons = defaultdict(list)


def isfunctional(csq):
    return any(c in ('stop_gained', 'stop_lost', 'start_lost', 'initiator_codon_variant', 'rare_amino_acid_variant',
                     'missense_variant', '5_prime_UTR_premature_start_codon_gain_variant', 'protein_altering_variant',
                     'splice_acceptor_variant', 'splice_donor_variant', 'frameshift_variant')
               for c in csq['Consequence'].split('&'))

# TODO: add a column that has all coverage values (% @10X for max(current.transcript.start, previous.end + 1):current.start
# TODO:                     ... same for GERP scores
print "#chrom\tstart\tend\taf\tfunctional\tgene\ttranscript\texon\timpact\tcnda_start\tcdna_end"
for i, v in enumerate(exac, start=1):
    if not (v.FILTER is None or v.FILTER == "PASS"):
        continue
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

    # skipping intronic
    if csq['Feature'] == '' or csq['EXON'] == '': continue

    cdna_start = csq['cDNA_position']
    if "-" in cdna_start:
        cdna_start, cdna_end = cdna_start.split("-")
    else:
        cdna_end = str(int(cdna_start) + len(v.REF))

    print "\t".join(map(str, (v.CHROM, v.start, v.end, af, int(isfunctional(csq)),
                              csq['SYMBOL'], csq['Feature'], csq['EXON'], csq['Consequence'],
                              cdna_start, cdna_end)))
