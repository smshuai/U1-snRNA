#!/usr/bin/env python
from pybedtools import BedTool
import pysam
import numpy as np
import os
import sys

## Arguments: miniBAM out_prefix
if (len(sys.argv)==3):
    bam_in = sys.argv[1]
    bed_path = sys.argv[2]
    out_prefix = 'out'
elif (len(sys.argv)==4):
    bam_in = sys.argv[1]
    bed_path = sys.argv[2]
    out_prefix = sys.argv[3]
else:
    sys.stderr.write("Usage: bam2fq_for_realign.py miniBAM U1_BED [OUT_PREFIX]")
    sys.exit(1)

mq_cut = 100  # reads < mapq_cut will be realigned

bam_tmp = os.path.splitext(bam_in)[0] + '.tmp.bam'
# Remove duplicates (1024)
pysam.view('-o', bam_tmp, '-uF', '1024', '-f', '1', bam_in, catch_stdout=False)
# Convert miniBAM to FASTQ
fq = pysam.fastq('-O', bam_tmp).split('\n')
# Convert FASTQ to dict
fq_dict = {}
for i in range(len(fq)):
    if i % 4 == 0:
        read_name = fq[i]
    elif i % 4 == 1:
        dna = fq[i]
    elif i % 4 == 3:
        qual = fq[i]
        # Record one read into dict
        fq_dict[read_name] = (dna, qual)
read_names = list(fq_dict.keys())
# Convert miniBAM to BED
U1U11 = BedTool(bed_path)
bam = BedTool(bam_tmp)
bed = bam.bam_to_bed().intersect(U1U11, wa=True, wb=True).to_dataframe()
bed['SELF_NAME'] = '@' + bed['name']
bed['QNAME'] = [i.split('/')[0] for i in bed['SELF_NAME']]  # QNAME as in BAM
bed['ORDER'] = [i.split('/')[1] for i in bed['SELF_NAME']]  # First or second in pair
# Add read names for mates
bed.loc[bed['ORDER']=='1', 'MATE_NAME'] = bed.loc[bed['ORDER']=='1', 'QNAME'] + '/2'
bed.loc[bed['ORDER']=='2', 'MATE_NAME'] = bed.loc[bed['ORDER']=='2', 'QNAME'] + '/1'
# Add some True/False indicators
bed['SELF'] = bed['SELF_NAME'].isin(read_names)  # is self in FASTQ?
bed['MATE'] = bed['MATE_NAME'].isin(read_names)  # is mate in FASTA?
bed['PAIRED'] = np.logical_and(bed['SELF'], bed['MATE'])  # Whether the pair can be found in the FASTQ
# Remove reads that are not in FASTQ
bed = bed[bed.SELF]
# Record names for re-align reads
keep = np.logical_and(bed['PAIRED'], bed['score'] < mq_cut)  # Realign reads with MAPQ < cut
use_qname = bed.loc[keep, ].QNAME.unique()
unpair_lowmq_qname = bed.loc[np.logical_and(np.logical_not(bed['PAIRED']), bed['score'] < mq_cut), ].SELF_NAME.unique()
# Print summary stat
print('Total number of reads:', len(read_names))
print('Total number of reads in BED:', bed['SELF'].sum())
print('Total number of paired reads in BED:', bed['PAIRED'].sum())
print('---> MAPQ>{}:'.format(mq_cut), np.logical_and(bed['PAIRED'], bed['score']>=mq_cut).sum(),
'and MAPQ<{}:'.format(mq_cut), np.sum(keep))
print('Total number of unpaired reads in BED:', np.logical_not(bed['PAIRED']).sum())
print('---> MAPQ>{}:'.format(mq_cut), np.logical_and(np.logical_not(bed['PAIRED']), bed['score']>=mq_cut).sum(),
      'and MAPQ<{}:'.format(mq_cut), np.logical_and(np.logical_not(bed['PAIRED']), bed['score']<mq_cut).sum())
print('Total Number of paired reads to be re-aligned:', len(use_qname)*2)
print('Total Number of unpaired reads to be re-aligned:', len(unpair_lowmq_qname))
# Make FASTQ for realignments
with open(out_prefix + '_1.fq', 'w') as f1:
    with open(out_prefix + '_2.fq', 'w') as f2:
        for qname in use_qname:
            f1.writelines("\n".join([qname, fq_dict[qname + '/1'][0], '+', fq_dict[qname + '/1'][1] + '\n']))
            f2.writelines("\n".join([qname, fq_dict[qname + '/2'][0], '+', fq_dict[qname + '/2'][1] + '\n']))
with open(out_prefix + '_unpair.fq', 'w') as f:
    for qname in unpair_lowmq_qname:
        f.writelines("\n".join([qname, fq_dict[qname][0], '+', fq_dict[qname][1] + '\n']))
# Clean up
os.remove(bam_tmp)
print('Done!')
