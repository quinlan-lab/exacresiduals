Regions of constraint are always defined by the bases between functional variants.

### First Step
1. use original, (not decomposed) VCF and take most deleterious allele.
2. union all exons for all transcripts in a chrom, gene pair (not UTRs). (analysis is no-longer transcript-based)
3. split at PASS, functional variants

### Second Pass
1. split at low-coverage regions inside exons (40 bases or entire exon)
2. split at self-chains

### Final Step
+ Add annotations (coverage, gerp, CpG) to each region.
+ Maybe add number of number of variants in non-PASS tranch / impact (func. v. non-func.)

If we implement like this, then we can add whatever we want to the
things that cause splits and it won't change the logic much because
we'll do a merge of anything that requires a split.

Most of the small pieces of the program logic will be unchanged:

+ read_coverage
+ is_functional
+ get_ranges
+ split_ranges

and these do the harder work. We'll just re-organize the procedural code
to avoid spiraling complexity.
