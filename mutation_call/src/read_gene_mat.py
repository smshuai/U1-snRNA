""" Male read x gene matrix for weight calculation
"""
#!/usr/bin/env python
# coding: utf-8
import pysam, sys, os
import numpy as np
import pandas as pd
from pybedtools import BedTool

# Input
if len(sys.argv) == 4:
    bam = sys.argv[1] # bowtie2 re-aligned & filtered bam (rbam for short)
    out_prefix = sys.argv[2]
    U1U11 = sys.argv[3]  # U1U11 w/ pseudo genes
else:
    sys.stderr.write('Incorrect arguments. Usage: postprocess_realign.py realign_BAM out_prefix U1_bed')
    sys.exit(1)
rbam_BT = BedTool(bam)  # rbam in BedTool
rbed_BT = rbam_BT.bam_to_bed(ed=True)
rbed_DF = rbed_BT.to_dataframe()
# Intersect with U1U11
rbed_U1U11 = rbed_BT.intersect(U1U11, wa=True, wb=True).to_dataframe(names=('chrom', 'start', 'end', 'name', 'score', 'strand', 'chrom2', 'start2', 'end2', 'gene', 'score2', 'strand2'))
# Make read x gene matrix
multimap_ct = rbed_DF.name.value_counts()  # number of times that each read aligns
ct_mat = rbed_U1U11.pivot_table(index='name', columns='gene', values='start', aggfunc=len)
ct_mat = ct_mat.fillna(0).astype(np.int)
ct_mat['num_align'] = multimap_ct[ct_mat.index]
ct_mat['U1U11_map'] = ct_mat.loc[:, ('RNU1-1', 'RNU1-2', 'RNU1-27P', 'RNU1-28P', 'RNU1-3', 'RNU1-4', 'RNVU1-18', 'RNU11')].sum(1)
ct_mat['Other'] = np.where(ct_mat.num_align>ct_mat.U1U11_map, 1, 0)
ct_mat.to_csv(out_prefix + '.mstat.csv')