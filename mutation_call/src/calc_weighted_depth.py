#!/usr/bin/env python
''' Calculate weighted depth based on re-alignments
Usage: calc_weighted_depth.py bam map_stat umapper U1_bed id out
'''
import pysam, re, sys, os
import pandas as pd
import numpy as np
import subprocess

ReIndel = re.compile('([\+\-])([0-9]+)([ACGTNacgtn]+)')
ReStart = re.compile('\^.')
ReEnd = re.compile('\$')

def process_args():
    if len(sys.argv) == 6:
        bam = sys.argv[1]
        U1U11 = sys.argv[2]
        donor_id = sys.argv[3]
        out_path = sys.argv[4]
        mstat = pd.read_csv(sys.argv[5], index_col='name')
        use_weight = 'yes'
        return bam, mstat, U1U11, donor_id, out_path, use_weight
    elif len(sys.argv) == 5:
        bam = sys.argv[1]
        U1U11 = sys.argv[2]
        donor_id = sys.argv[3]
        out_path = sys.argv[4]
        use_weight = 'no'
        return bam, None, U1U11, donor_id, out_path, use_weight
    else:
        sys.stderr.write('Usage: calc_weighted_depth.py BAM U1_bed ID OUT [map_stat]\n')
        sys.stderr.write('\t- BAM: miniBAM file\n')
        sys.stderr.write('\t- U1_bed: BED file for U1/U11 genes\n')
        sys.stderr.write('\t- ID: sample identifier, also used as output prefix\n')
        sys.stderr.write('\t- OUT: path for output file\n')
        sys.stderr.write('\t- map_stat: (optional) Mapping statistics from post-processing step\n')
        sys.exit(1)

def get_weights(depth, readBar, ct):
    ''' Use readBar (array of read names), ct (number of times a read maps) to get weights for each read
    '''
    if depth == 0:
        return
    weights = []
    for read in readBar:
        if read in ct.index:
            w = 0 if ct.loc[read, 'Other'] > 0 else ct.loc[read, 'U1U11_map']  # reads mapped to other regions are ignored
            weights.append(int(w))
        else:
            sys.stderr('Read not in umap or mmap:', read)
            sys.exit(1)
    assert len(weights) == len(readBar)
    return(weights)

def prepWeightedCountCheck(baseBar, weights):
    insertion_vb_p, insertion_vb_n, deletion_vb_p, deletion_vb_n = 0, 0, 0, 0
    # fix indels in baseBar
    deleted = 0
    iter = ReIndel.finditer(baseBar)
    for m in iter:
        site = m.start()
        indelType = m.group(1)
        indelSize = m.group(2)
        varChar = m.group(3)[0:int(indelSize)]
        strand = '+' if varChar.isupper() else '-'
        if indelType == "+":
            if strand == "+":
                insertion_vb_p += 1
            else:
                insertion_vb_n += 1
        else:
            if strand == "+":
                deletion_vb_p += 1
            else:
                deletion_vb_n += 1
        baseBar = baseBar[0:(site - deleted)] + baseBar[(site + int(indelSize) + len(indelSize) + 1 - deleted):]
        deleted += 1 + len(indelSize) + int(indelSize)
    # Remove start and end info in baseBar
    baseBar = ReStart.sub('', baseBar)
    baseBar = ReEnd.sub('', baseBar)
    # error check
    assert len(baseBar) == len(weights), "lengths of bases and weights are different!"
    base_num_weighted = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "a": 0, "c": 0, "g": 0, "t": 0, "n": 0}
    base_num = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "a": 0, "c": 0, "g": 0, "t": 0, "n": 0}
    raw_dict = {"A": [], "C": [], "G": [], "T": [], "N": [], "a": [], "c": [], "g": [], "t": [], "n": []}
    # rely on -Q option in samtools; base quality check is removed
    for base, w in zip(baseBar, weights):
        if base in "ATGCNatgcn":
            base_num_weighted[base] += 1/w if w>0 else 0  # the most importane line; weight each base by its times of mapping
            base_num[base] += 1 if w>0 else 0  # unweighted version
            raw_dict[base].append(str(w))
    raw_dict = pd.Series(raw_dict).apply(lambda x: ",".join(x))
    return base_num_weighted, base_num, raw_dict, insertion_vb_n, insertion_vb_p, deletion_vb_n, deletion_vb_p


def varWeightedCountCheck(var, base_num_weighted, base_num, weighted):
    if weighted:
        # Depth is calcuated with weighted version
        depth_p = base_num_weighted["A"] + base_num_weighted["C"] + base_num_weighted["G"] + base_num_weighted["T"] + base_num_weighted["N"]
        depth_n = base_num_weighted["a"] + base_num_weighted["c"] + base_num_weighted["g"] + base_num_weighted["t"] + base_num_weighted["n"]
    else:
        # unweighted
        depth_p = base_num["A"] + base_num["C"] + base_num["G"] + base_num["T"] + base_num["N"]
        depth_n = base_num["a"] + base_num["c"] + base_num["g"] + base_num["t"] + base_num["n"]
    misMatch_p = 0
    misMatch_n = 0
    # Mismatch depth is calcuated with unweighted version
    if var in "ACGTacgt":
        # variant is SNV
        misMatch_p = base_num[var.upper()]
        misMatch_n = base_num[var.lower()]
    elif var[0] == "+":
        # variant is INS
        misMatch_p, misMatch_n = insertion_vb_p, insertion_vb_n
    elif var[0] == "-":
        # variant is DEL
        misMatch_p, misMatch_n = deletion_vb_p, deletion_vb_n
    else:
        sys.stderr.write(var + ": input var has wrong format!")
        sys.exit(1)
    return [int(misMatch_p), int(np.ceil(depth_p)), int(misMatch_n), int(np.ceil(depth_n))]


def run_subprocess(cmd):
    exit_status = subprocess.call(cmd, shell=True)
    if exit_status is 1:
        sys.stderr.write('Failed:', cmd)
        sys.exit(1)


def pileup2dict(pileup, mate_order, results, mstat):
    with open(pileup) as f:
        for pp in f:
            (chrom, pos, _, depth, baseBar, _, readBar) = pp.strip().split('\t')
            depth = int(depth)
            if depth > 0:
                readBar = [i + mate_order for i in readBar.split(',')]
                weights = get_weights(depth, readBar, mstat)
                base_num_weighted, base_num, raw_dict, insertion_vb_n, insertion_vb_p, deletion_vb_n, deletion_vb_p = prepWeightedCountCheck(baseBar, weights)
                for var in 'ACGT':
                    key = (chrom, pos, var)
                    counts = varWeightedCountCheck(var, base_num_weighted, base_num, weighted=True)
                    # use multiple pileup files, so need to perform addition
                    results[key] = np.add(results.get(key, [0, 0, 0, 0]), counts).tolist()
            else:
                for var in 'ACGT':
                    key = (chrom, pos, var)
                    counts = [0, 0, 0, 0]
                    results[key] = np.add(results.get(key, [0, 0, 0, 0]), counts).tolist()
    return results


def pileup2dict_noweight(pileup, results):
    with open(pileup) as f:
        for pp in f:
            (chrom, pos, _, depth, baseBar, qualBar) = pp.strip().split('\t')
            depth = int(depth)
            if depth > 0:
                weights = [1 for i in qualBar]
                base_num_weighted, base_num, raw_dict, _, _, _, _ = prepWeightedCountCheck(baseBar, weights)
                for var in 'ACGT':
                    key = (chrom, pos, var)
                    counts = varWeightedCountCheck(var, base_num_weighted, base_num, weighted=False)
                    results[key] = np.add(results.get(key, [0, 0, 0, 0]), counts).tolist()
            else:
                for var in 'ACGT':
                    key = (chrom, pos, var)
                    counts = [0, 0, 0, 0]
                    results[key] = np.add(results.get(key, [0, 0, 0, 0]), counts).tolist()
    return results


if __name__ == '__main__':
    # Get arguments
    bam, mstat, U1U11, donor_id, out_path, use_weight = process_args()
    results = dict()
    if use_weight is 'yes':
        # Make two pileup files (/1, /2)
        cmd1 = 'samtools mpileup -l {} -d 1000000 --rf 64 --ff UNMAP,QCFAIL -a --output-QNAME -Q 15 {} > {}'.format(U1U11, bam, bam + '_1.pileup')
        cmd2 = 'samtools mpileup -l {} -d 1000000 --rf 128 --ff UNMAP,QCFAIL -a --output-QNAME -Q 15 {} > {}'.format(U1U11, bam, bam + '_2.pileup')
        run_subprocess(cmd1)
        run_subprocess(cmd2)
        # Process /1 pileup
        results = pileup2dict(bam + '_1.pileup', '/1', results, mstat)
        # Process /2 pileup
        results = pileup2dict(bam + '_2.pileup', '/2', results, mstat)
    else:
        # No weight
        cmd = 'samtools mpileup -l {} -d 1000000 --ff UNMAP,QCFAIL -a -Q 15 {} > {}'.format(U1U11, bam, bam + '.pileup')
        run_subprocess(cmd)
        # Process pileup
        results = pileup2dict_noweight(bam + '.pileup', results)
    res_DF = pd.DataFrame(results).T.reset_index()
    res_DF.columns = ('chrom', 'pos', 'alt', 'varP', 'dpP', 'varN', 'dpN')
    res_DF['id'] = donor_id
    # Save result
    res_DF.to_csv(out_path, sep='\t', index=False)
    # Clean up
    if use_weight is 'yes':
        os.remove(bam + '_1.pileup')
        os.remove(bam + '_2.pileup')
    else:
        os.remove(bam + '.pileup')
