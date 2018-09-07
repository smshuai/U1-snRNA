#!/usr/bin/env python
# coding: utf-8
import pysam, sys, os
import numpy as np
from pybedtools import BedTool

mq_cut = 100

# Input
if len(sys.argv) == 5:
    realnBAM = sys.argv[1] # bowtie2 re-aligned bam (rbam for short)
    oriBAM = sys.argv[2]  # original miniBAM (obam)
    out_prefix = sys.argv[3]
    U1U11 = sys.argv[4]  # U1U11 w/ pseudo genes
else:
    sys.stderr.write('Incorrect arguments. Usage: postprocess_realign.py realign_BAM original_BAM out_prefix U1_bed')
    sys.exit(1)
####
# Remove non-optimal multiple alignments
####
rbam_BT = BedTool(realnBAM)  # rbam in BedTool
rbed_DF = rbam_BT.bam_to_bed(tag='AS').to_dataframe()
# Bowtie2 end-to-end returns negative scores so bedtools add 256 to the AS for negative values.
# Therefore, 0 score should actually be fixed to 256
rbed_DF.loc[rbed_DF.score==0, 'score'] = 256
# Only keep multi-alignments with max AS per read
max_score = rbed_DF.pivot_table(index='name', values='score', aggfunc=max).reset_index()
rbed_DF = rbed_DF.set_index(['name', 'score']).sort_index()
keep = max_score.set_index(['name', 'score']).index.values.tolist()
rbed_DF = rbed_DF.loc[keep, ].reset_index()
rbed_DF = rbed_DF.loc[:,('chrom', 'start', 'end', 'name', 'score', 'strand')]
multimap_ct = rbed_DF.name.value_counts()  # number of times that each read aligns
rbed_BT = BedTool.from_dataframe(rbed_DF)  # create filtered version of rbed
# Intersect with core U1U11
rbed_BT = rbed_BT.intersect(U1U11, wa=True, wb=True)
rbed_DF = rbed_BT.to_dataframe(names=('chrom', 'start', 'end', 'name', 'score', 'strand', 'chrom2', 'start2', 'end2', 'gene', 'score2', 'strand2'))
# Make read x gene matrix
ct_mat = rbed_DF.pivot_table(index='name', columns='gene', values='strand', aggfunc=len)
ct_mat['num_align'] = multimap_ct[ct_mat.index]
ct_mat.fillna(0, inplace=True)
ct_mat['U1_map'] = ct_mat.loc[:, ('RNU1-1', 'RNU1-2', 'RNU1-27P', 'RNU1-28P', 'RNU1-3', 'RNU1-4', 'RNVU1-18')].sum(1)
ct_mat['Other'] = np.where(ct_mat.num_align>ct_mat.U1_map, 1, 0)
####
# Make SAM for EB call
####
# Subset bam first to save some time
rbam_U1U11 = rbam_BT.intersect(rbed_BT, r=True, f=1)
obam_U1U11 = BedTool(oriBAM).intersect(U1U11)
umap_reads = []
with open(out_prefix + '.processed.sam', 'w') as out:
    # Add original unique mapping reads to OUT SAM
    with pysam.Samfile(obam_U1U11.fn, 'rb') as obam:
        out.write(str(obam.header))
        for read in obam:
            if read.mapq >= mq_cut:
                read_name = read.qname + '/2' if read.is_read2 else read.qname + '/1'
                umap_reads.append(read_name)  # store names of unique mapping reads
                out.write(read.to_string() + '\n')
    # Remove unique mapping reads from rbed
    rbed_DF = rbed_DF.loc[np.logical_not(rbed_DF.name.isin(umap_reads))]
    # Add re-aligned multiple mapping reads to OUT SAM
    with pysam.Samfile(rbam_U1U11.fn, 'rb') as rbam:
        for read in rbam:
            read_name = read.qname + '/2' if read.is_read2 else read.qname + '/1'
            if (not read.is_read2) and (not read.is_read1):
                read_name = read.qname  # unpaired reads
            if read_name in rbed_DF.name.values:
                starts = rbed_DF.loc[rbed_DF.name==read_name, 'start'].values
                if read.reference_start in starts:
                    out.write(read.to_string() + '\n')
# Sort SAM
pysam.sort("-o", out_prefix + '_processed.bam', out_prefix + '.processed.sam', catch_stdout=False)
# Remove umappers
ct_mat = ct_mat.loc[np.logical_not(ct_mat.index.isin(umap_reads))]
ct_mat.to_csv(out_prefix + '_mstat.csv')
# rbed_DF.to_csv('./test_rbed.csv', index=False)
with open(out_prefix + '_umap.txt', 'w') as f:
    for read in umap_reads:
        f.write(read + '\n')
# Clean up
os.remove(out_prefix + '.processed.sam')
print('Use {} unique mapping reads and {} multiple mapping reads'.format(len(umap_reads), rbed_DF.shape[0]))
