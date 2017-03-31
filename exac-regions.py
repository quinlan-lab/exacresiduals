from __future__ import print_function

# ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/ExAC.r0.3.sites.vep.vcf.gz
VCF_PATH = "data/ExAC.r0.3.sites.vt.vep.vcf.gz" #"toyexac.vcf.gz" #"data/gnomad.exomes.r2.0.1.sites.vcf.gz"

# ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
GTF_PATH = "data/Homo_sapiens.GRCh37.75.gtf.gz" #"toyexons.gtf.gz" #"data/Homo_sapiens.GRCh37.75.gtf.gz"

# ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/coverage
COVERAGE_PATH = "data/"

# from UCSC. see data/get-chain.py, pipe output to sort -k1,1 -k2,2n | uniq | bgzip -c > data/self-chains.gt90.bed.gz
SELF_CHAINS = "data/self-chains.gt90.bed.gz"

# from UCSC.  data/segmental.bed.gz

SEGDUPS = "data/segmental.bed.gz"

FASTA_PATH = "/uufs/chpc.utah.edu/common/home/u6000771/Data/data/hs37d5.fa"


import sys
import itertools as it
import operator

import numpy as np
import utils as u
from cyvcf2 import VCF
from pyfaidx import Fasta

import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-n", "--nosingletons", help="if you do NOT want singletons", action="store_true", default=False)
parser.add_argument("-v", "--varflag", help="if you want separation by variant flags", action="store_true", default=False)
args=parser.parse_args()
nosingletons=args.nosingletons
varflag=args.varflag

zip = it.izip


exac = VCF(VCF_PATH)
kcsq = exac["CSQ"]["Description"].split(":")[1].strip(' "').split("|")

#exac = exac("2:112538945-112551053")


fasta = Fasta(FASTA_PATH, read_ahead=10000, as_raw=True)

header = "chrom\tstart\tend\taf\tfunctional\tgene\ttranscript\texon\timpact\tvstart\tvend\tn_bases\tcg_content\tranges\tcoverage\tposns\tvarflag"
print("#" + header)
keys = header.split("\t")
global mranges, splitter


def merge_rows(rows):
    """
    >>> merge_rows([dict(gene='ABC', vstart=1, vend=3), dict(gene='ABC', vstart=2, vend=5), dict(gene='ABC', vstart=4, vend=6)])
    [{'vstart': 1, 'vend': 6, 'gene': 'ABC'}]
    >>> merge_rows([dict(gene='ABC', vstart=1, vend=6), dict(gene='ABC', vstart=2, vend=6), dict(gene='ABC', vstart=4, vend=6)])
    [{'vstart': 1, 'vend': 6, 'gene': 'ABC'}]
    >>> merge_rows([dict(gene='ABC', vstart=1, vend=2), dict(gene='ABC', vstart=1, vend=4), dict(gene='ABC', vstart=1, vend=6)])
    [{'vstart': 1, 'vend': 6, 'gene': 'ABC'}]
    >>> merge_rows([dict(gene='ABC', vstart=1, vend=2), dict(gene='Easy as 123', vstart=1, vend=4), dict(gene='Easy as 123', vstart=1, vend=6)])
    [{'vstart': 1, 'vend': 2, 'gene': 'ABC'}, {'vstart': 1, 'vend': 6, 'gene': 'Easy as 123'}]
    >>> merge_rows([dict(gene='ABC', vstart=1, vend=2), dict(gene='ABC', vstart=2, vend=4), dict(gene='ABC', vstart=4, vend=6)])
    [{'vstart': 1, 'vend': 2, 'gene': 'ABC'}, {'vstart': 2, 'vend': 4, 'gene': 'ABC'}, {'vstart': 4, 'vend': 6, 'gene': 'ABC'}]
    >>> merge_rows([dict(gene='ABC', vstart=1, vend=2), dict(gene='Easy as', vstart=1, vend=4), dict(gene='123', vstart=1, vend=6)])
    [{'vstart': 1, 'vend': 2, 'gene': 'ABC'}, {'vstart': 1, 'vend': 4, 'gene': 'Easy as'}, {'vstart': 1, 'vend': 6, 'gene': '123'}]
    >>> merge_rows([dict(gene='ABC', vstart=1, vend=2), dict(gene='ABC', vstart=1, vend=6), dict(gene='ABC', vstart=1, vend=4)])
    [{'vstart': 1, 'vend': 6, 'gene': 'ABC'}]
    """
    new_rows = [rows[0]]
    for row in rows[1:]:
        if new_rows[-1]['gene'] == row['gene'] and row['vstart'] < new_rows[-1]['vend'] and row['vend'] > new_rows[-1]['vend']:
            new_rows[-1]['vend'] = row['vend']
        elif new_rows[-1]['gene'] != row['gene'] or new_rows[-1]['vend'] <= row['vstart']:
            new_rows.append(row)
    return new_rows

def separate_ranges(ranges, varflags): # for putting each VARTRUE range in its own range list, so that coverage and "cg_content" are only for true regions
    """
    >>> separate_ranges([[(26782, 26890)], [(45349, 45487), (58320, 58416), (58687, 58777), (60611, 60618)]], [['VARFALSE'], ['VARFALSE', 'VARFALSE', 'VARFALSE', 'VARTRUE']])
    ([[(26782, 26890)], [(45349, 45487), (58320, 58416), (58687, 58777)], [(60611, 60618)]], [['VARFALSE'], ['VARFALSE', 'VARFALSE', 'VARFALSE'], ['VARTRUE']])
    >>> separate_ranges([[(26782, 26890)], [(45349, 45487), (58320, 58416), (58687, 58777), (60611, 60618), (60422, 60443)]], [['VARFALSE'], ['VARFALSE', 'VARFALSE', 'VARFALSE', 'VARTRUE', 'VARTRUE']])
    ([[(26782, 26890)], [(45349, 45487), (58320, 58416), (58687, 58777)], [(60611, 60618), (60422, 60443)]], [['VARFALSE'], ['VARFALSE', 'VARFALSE', 'VARFALSE'], ['VARTRUE', 'VARTRUE']])
    """
    newranges=[]; newvarflags=[]
    for rangelist, flaglist in zip(ranges, varflags):
        fr, fv, tr, tv = [], [], [], []
        for r, v in zip(rangelist, flaglist):
            if v=='VARFALSE':
                fr.append(r); fv.append(v) 
            if v=='VARTRUE':
                tr.append(r); tv.append(v)
        if fr and fv:
            newranges.append(fr)
            newvarflags.append(fv)
        if tr and tv:
            newranges.append(tr)
            newvarflags.append(tv)
  
    return newranges, newvarflags

import doctest
res = doctest.testmod(verbose=0)
if res.failed != 0:
    sys.exit(1)

for chrom, viter in it.groupby(exac, operator.attrgetter("CHROM")):
    rows = []
    print("reading chrom", file=sys.stderr)

    fa = fasta[chrom]
    coverage_array = u.read_coverage(chrom, length=len(fa), cov=10,
            path=COVERAGE_PATH)

    gene_exon_starts, gene_exon_ends, splitters = u.read_exons("|tabix {gtf} {chrom}"
                                                            .format(chrom=chrom,
                                                                gtf=GTF_PATH),
                                                            chrom,
                                                            coverage_array,
                                            "|tabix {bed} {chrom}".format(chrom=chrom, bed=SELF_CHAINS),"|tabix {bed} {chrom}".format(chrom=chrom, bed=SEGDUPS))

    print(chrom, file=sys.stderr)
    for v in viter:
        if not (v.FILTER is None or v.FILTER == "PASS"):
            continue
        info = v.INFO
        try:
            csqs = [dict(zip(kcsq, c.split("|"))) for c in info['CSQ'].split(",")]
        except KeyError:
            continue
        # NOTE: using max here for alternates to be conservative
        ac = info['AC_Adj']
        if not isinstance(ac, (int, long)):
            ac = max(ac)
        af = ac / float(info['AN_Adj'] or 1)
        if ac == 1: #self-explanatory, but filters out singletons
            if nosingletons: continue
        # NOTE: not requiring canonical or requiring the csq to match the
        # particular alt that we chose.
        for csq in (c for c in csqs if c['BIOTYPE'] == 'protein_coding'): # getting duplicate rows because of this, wastes memory and potentially compute time, could remove and replace with just if isfunctional, add to rows then move on?
            # skipping intronic
            if csq['Feature'] == '' or csq['EXON'] == '' : continue #or csq['cDNA_position'] == '': continue
            if not u.isfunctional(csq): continue

            if csq['cDNA_position']:
                cdna_start, cdna_end = u.get_cdna_start_end(csq['cDNA_position'], v)
            else:
                #print(v,file=sys.stderr)
                cdna_start, cdna_end = 'na','na' # apparently, sometimes ENSEMBL doesn't annotate splice_donor_variant&coding_sequence_variant combinations with cdna coords

            rows.append(dict(chrom=v.CHROM, vstart=v.start, vend=v.end, af=af,
                functional=int(u.isfunctional(csq)),
                gene=csq['SYMBOL'], transcript=csq['Feature'], exon=csq['EXON'],
                impact=csq['Consequence'],
                cdna_start=cdna_start,   cdna_end=cdna_end))

        

    # now we need to sort and then group by gene so we know the gaps.
    rows.sort(key=operator.itemgetter('gene', 'vstart', 'vend'))

#TODO:we're only using vstart. maybe we should be using end, normalizing and decomposing?
# we may be misrepresenting deletions as in the case of vstart, vend, pos: 145838622 145838695 145838623; shouldn't this whole area not be covered?
    
    rows = merge_rows(rows)

    out = []
    for chrom_gene, trows in it.groupby(rows, lambda row: (row['chrom'], row['gene'])):
        exon_starts = gene_exon_starts[chrom_gene]
        exon_ends = gene_exon_ends[chrom_gene]
        last = exon_starts[0]
        splitter = splitters.get(chrom_gene, None)

        for i, row in enumerate(trows, start=1):
            # istart and iend determine if we need to span exons.

            assert row['vstart'] <= exon_ends[-1], (row, exon_ends) # maybe use POS instead of vstart, so we can normalize and decompose?; should i check if end is less?
            row['vstart']=row['vstart']+1 # vstart is bed format variant coordinate, still true maybe use POS instead of vstart?
            last2 = last
            if varflag:
                mranges, last, varflags = u.get_ranges(last, row['vstart'], row['vend'], exon_starts, exon_ends, row['chrom'])#TODO: fix get_ranges to do what split_ranges does, and land behind vend because it ends at vstart
 #           print (last, row['vstart'], row['vend'], mranges, splitter, varflags)
            else:
                 mranges, last, varflags = u.get_ranges_w_variant(last, row['vstart'], row['vend'], exon_starts, exon_ends, row['chrom'])
            mranges2, varflags2 = u.split_ranges(mranges, splitter, varflags)
            if varflag:
                mranges2, varflags2 = separate_ranges(mranges2, varflags2)
            for ranges, vf in zip(mranges2, varflags2):
                row['coverage'] = ",".join(",".join(u.floatfmt(g) for g in coverage_array[s:e]) for s, e in ranges)
                row['posns'] = list(it.chain.from_iterable([range(s+1, e+1) for s, e in ranges])) # since range is not inclusive at the end add +1, need to add +1 to start
                row['ranges'] = ["%d-%d" % (s, e) for s, e in ranges]
                row['varflag'] = ",".join(vf)
                #print (row)
                #print (last,row['vstart'], "last n' vstart")
                #print (ranges, row['ranges'])
                #print (exon_starts, exon_ends, "starts n' ends")
                #print (splitter, "splitta!")
                seqs = [fa[s:e] for s, e in ranges]
                # this can happen for UTR variants since we can't really get
                # anything upstream of them.
                if row['posns'] == []:  # UTR:
                    p = row['vstart']
                    row['coverage'] = ",".join(u.floatfmt(g) for g in coverage_array[p:p+1])
                    row['posns'] = [p]

                # post-hoc sanity check
                exon_bases = set(it.chain.from_iterable(range(s, e) for s, e in zip(exon_starts, exon_ends)))
                ranges = set(it.chain.from_iterable(range(int(x[0]), int(x[1])) for x in (z.split("-") for z in row['ranges'])))
                m = len(ranges - exon_bases)
                if m > len(row['ranges']):
                    #print (last2, row['vstart'], row['vend'], exon_starts, exon_ends)
                    #print (row['ranges'], ranges - exon_bases)
                    print(last, row['vstart'], row['ranges'], len(ranges -
                        exon_bases), zip(exon_starts, exon_ends),
                        file=sys.stderr)

                # start or end? if we use end then can have - diff.
                row['ranges'] = ",".join(row['ranges'])
                row['n_bases'] = len(row['posns'])
                row['start'] = str(min(row['posns'])-1) #I put -1 because I am not including the position of the start coordinate, as it is 0-based.  however, I want to make it the start coordinate
                row['end'] = str(max(row['posns'])) # base on ranges
                row['posns'] = ",".join(map(str, row['posns']))
                row['cg_content'] = u.floatfmt(np.mean([u.cg_content(s) for s in seqs]))
                if row['cg_content'] == 'nan':
                    row['cg_content'] = '0'
                # we are re-using the dict for each loop so force a copy.
                try:
                    if row['ranges']:
                        endrange=int(row['ranges'].split('-')[-1])
                        if last < endrange:
                            last = endrange  #so we can start at where the last range ended
                        out.append(dict(row))
                except KeyError:
                    pass
            #else:
            #    row['ranges']=row['start']+"-"+row['end']
            #    row['end']=str(int(row['end'])+1)
            #    out.append(dict(row))
            # last = row['vstart']

        for exon_ending in exon_ends: # when variant is the last variant in a gene, it fails to make regions for the end of the gene, thus this block of code
            try:
                if mranges[-1][-1] > exon_ending: # do if mranges[-1][-1] < exon_ending or mranges is empty; if it is equal, we want it because it could be a 1 bp region
                    continue
            except IndexError:
                pass 
#            print (last, exon_ends[-1], row['vend'], exon_starts, exon_ends, row['chrom'])
            mranges, last, varflags = u.get_ranges(last, exon_ends[-1]+1, exon_ends[-1]+1, exon_starts, exon_ends, row['chrom']) #TODO: fix vend?
            
            mranges2, varflags2 = u.split_ranges(mranges, splitter, varflags)
            if varflag:
                mranges2, varflags2 = separate_ranges(mranges2, varflags2)
            for ranges, vf in zip(mranges2, varflags2):
                row['coverage'] = ",".join(",".join(u.floatfmt(g) for g in coverage_array[s:e]) for s, e in ranges)
                row['posns'] = list(it.chain.from_iterable([range(s+1, e+1) for s, e in ranges])) #range is not inclusive at the end, need to add +1 to s
                row['ranges'] = ["%d-%d" % (s, e) for s, e in ranges]
                row['varflag'] = ",".join(vf)
                #print (row)
                #print (last,row['vstart'], "last n' vstart")
                #print (ranges, row['ranges'])
                #print (exon_starts, exon_ends, "starts n' ends")
                #print (splitter, "splitta!")
                seqs = [fa[s:e] for s, e in ranges]
                # this can happen for UTR variants since we can't really get
                # anything upstream of them.
                if row['posns'] == []:  # UTR:
                    p = row['vstart']
                    row['coverage'] = ",".join(u.floatfmt(g) for g in coverage_array[p:p+1])
                    row['posns'] = [p]
                # post-hoc sanity check
                exon_bases = set(it.chain.from_iterable(range(s, e) for s, e in zip(exon_starts, exon_ends)))
                ranges = set(it.chain.from_iterable(range(int(x[0]), int(x[1])) for x in (z.split("-") for z in row['ranges'])))
                m = len(ranges - exon_bases)
                if m > len(row['ranges']):
                    print(last, row['vstart'], row['ranges'], len(ranges -
                        exon_bases), zip(exon_starts, exon_ends),
                        file=sys.stderr)

                # start or end? if we use end then can have - diff.
                row['ranges'] = ",".join(row['ranges'])
                row['n_bases'] = len(row['posns'])
                row['start'] = str(min(row['posns'])-1) #I put -1 because I am not including the position of the start coordinate, as it is 0-based.  however, I want to make it he start coordinate 
                row['end'] = str(max(row['posns']))
                row['posns'] = ",".join(map(str, row['posns']))
                row['cg_content'] = u.floatfmt(np.mean([u.cg_content(s) for s in seqs]))
                if row['cg_content'] == 'nan':
                    row['cg_content'] = '0'
                # we are re-using the dict for each loop so force a copy.
                try:
                    if row['ranges']:
                        last = int(row['ranges'].split('-')[-1]) #so we can start at where the last range ended
                        out.append(dict(row))
                except KeyError:
                    pass 

    # still print in sorted order
    out.sort(key=operator.itemgetter('start'))
    last = (None, None)
    for d in out:
        key = d['start'], d['end'], d['gene'] # added d['gene'] so it doesn't throw away longer regions without variants in a different gene overlapping the same genome space
        if key == last:
            continue
        last = key
        print("\t".join(map(str, (d[k] for k in keys))))
