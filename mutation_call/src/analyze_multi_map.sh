#!/usr/bin/env bash
# Author: Shimin Shuai
# Date: 2018-09-06
# Contact: shimin.shuai@mail.utoronto.ca
###

### Input arguments
miniBAM=$1  # 'miniBAM/e6cdf28912d360a2b947f5911c79e214.gencode.v19.U1U11.flank1k.nochr.bed.bam'
donorID=$2  # donor/patient identifier
OUT=$3  # output prefix

### Required reference files; might need modification
BW2IX="/u/sshuai/sshuai/ref/bowtie2index/genome"  # set your bowtie2 index
U1BED="../ref/gencode.v19.U1U11.nochr.bed"  # U1 bed with all variants and pseudogenes
U1CORE="../ref/gencode.v19.U1U11.nochr.nopseduo.bed"  # Core U1 genes (7 U1 + 1 U11) 
# Software and references
BW2=bowtie2  # bowtie2
src_dir="./"  # Directory for scripts if not in PATH
BAM2FQ="$src_dir/bam2fq_for_realign.py"  # script for converting bam to fastq
POSTPROCESS="$src_dir/filter_realn.py"  # script for postprocessing
READGENE="$src_dir/read_gene_mat.py"
CALCDP="$src_dir/calc_weighted_depth.py"  # script used to calcuate weighted depth


#################
### PART 1
### Realign with bowtie2
#################
echo "> Make FASTQ for realignments"
python $BAM2FQ $miniBAM $U1BED $OUT
echo "> Run Bowtie2"
$BW2 --score-min L,-0.3,-0.3 --no-mixed --no-discordant -k 100 --very-sensitive -p 4 -x $BW2IX -1 ${OUT}_1.fq -2 ${OUT}_2.fq -S ${OUT}.realn.sam
samtools view -bh ${OUT}.realn.sam > ${OUT}.realn.bam
#################
### PART 2
### Filter re-aligned BAM
#################
echo "> Filter BAM"
python $POSTPROCESS ${OUT}.realn.bam $OUT
echo "> Make read X gene matrix"
python $READGENE ${OUT}.filtered.realn.bam $OUT $U1BED
#################
### PART 3
### Make weighted depth table
#################
echo "> Make unweighted depth table"
python $CALCDP ${OUT}.filtered.realn.bam  $U1CORE $donorID ${OUT}.unweighted.dp.tsv
echo "> Make weighted depth table"
python $CALCDP ${OUT}.filtered.realn.bam $U1CORE $donorID ${OUT}.weighted.dp.tsv ${OUT}.mstat.csv
## Clean up
rm ${OUT}*.fq ${OUT}.realn.sam

