This directory contains scripts used to call somatic mutations in U1 snRNA.

## STEP 0: Make miniBAMs from BWA-MEM alignments
Note that BWA-MEM must be run with `-t 0` to record all low MQ reads. Then for each BAM,
samtools can be used to extract reads mapped to U1 regions and make miniBAMs.
```bash
$ samtools view -L $BED -hu $BWA_BAM | samtools sort > $miniBAM
```
`$BED` is the BED file for U1 genes and pesudogenes and their flanking sequences. The human version (GENCODE v19; GRCh37) of `$BED` is provided.

## STEP 1: Analyze multiple mapping patterns
This step has three parts: 1. realign with Bowtie2; 2. Filter realigned BAM; 3. Make REF and VAR depth table.

A warpper script is provided, which takes three positional arguments:
1. one miniBAM from **STEP 0**; 2. the identifier of the donor; 3. the output BAM path
```bash
$ analyze_multi_map.sh $miniBAM $DONOR_ID $BW2_BAM
```
