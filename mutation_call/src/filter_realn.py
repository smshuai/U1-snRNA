#!/usr/bin/env python
# coding: utf-8
import pysam, sys, os
import numpy as np
import pandas as pd
from pybedtools import BedTool

# Input
if len(sys.argv) == 3:
    realnBAM = sys.argv[1] # bowtie2 re-aligned bam (rbam for short)
    out_prefix = sys.argv[2]
else:
    sys.stderr.write('Incorrect arguments. Usage: postprocess_realign.py realign_BAM out_prefix')
    sys.exit(1)
####
# Remove non-optimal multiple alignments
####
rbam_BT = BedTool(realnBAM)  # rbam in BedTool
cnames = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'ed', 'strand1', 'strand2']
rbed_DF = pd.read_table(rbam_BT.bam_to_bed(ed=True, mate1=True, bedpe=True).fn, header=None, names=cnames)
rbed_DF = rbed_DF.loc[rbed_DF['chrom1']!='.']
rbed_DF = rbed_DF.loc[rbed_DF['chrom2']!='.']
# Only keep multi-alignments with min ED per pair
min_ed = rbed_DF.pivot_table(index='name', values='ed', aggfunc=min).reset_index()
rbed_DF = rbed_DF.set_index(['name', 'ed']).sort_index()
keep = min_ed.set_index(['name', 'ed']).index.values.tolist()
rbed_DF = rbed_DF.loc[keep, ].reset_index()
# bedpe to BED
mate1 = rbed_DF.loc[:, ('chrom1', 'start1', 'end1')]
mate1['name'] = rbed_DF['name'] + '/1'
mate2 = rbed_DF.loc[:, ('chrom2', 'start2', 'end2')]
mate2['name'] = rbed_DF['name'] + '/2'
mate1.columns = ('chrom', 'start', 'end', 'name')
mate2.columns = ('chrom', 'start', 'end', 'name')
mates = mate1.append(mate2)
####
# Make SAM for EB call
####
with open(out_prefix + '.filtered.realn.sam', 'w') as out:
    # Write header to output sam
    out.write(pysam.view("-H", realnBAM))
    # Write filtered reads to OUT SAM
    with pysam.Samfile(realnBAM, 'rb') as rbam:
        for read in rbam:
            read_name = read.qname + '/2' if read.is_read2 else read.qname + '/1'
            if (not read.is_read2) and (not read.is_read1):
                read_name = read.qname  # unpaired reads
            if read_name in mates.name.values:
                starts = mates.loc[mates.name==read_name, 'start'].values
                if read.reference_start in starts:
                    out.write(read.to_string() + '\n')
# Sort SAM
pysam.sort("-o", out_prefix + '.filtered.realn.bam', out_prefix + '.filtered.realn.sam', catch_stdout=False)
# Clean up
os.remove(out_prefix + '.filtered.realn.sam')