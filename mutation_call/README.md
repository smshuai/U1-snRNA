This directory contains scripts used to call somatic mutations in U1 snRNA. The required inputs are normal-tumour paired WGS data.

## STEP 0: Make miniBAMs from BWA-MEM alignments
Note that BWA-MEM must be run with `-t 0` to record all low MQ reads. Then for each BAM,
samtools can be used to extract reads mapped to U1 regions and make miniBAMs.
```bash
$ samtools view -L $BED -hu $BWA_BAM | samtools sort > $miniBAM
```
`$BED` is the BED file for U1 genes and pesudogenes and their flanking sequences. The human version (GENCODE v19; GRCh37) of `$BED` is provided ([./ref/gencode.v19.U1U11.nochr.flank1k.bed](./ref/gencode.v19.U1U11.nochr.flank1k.bed)).

## STEP 1: Analyze multiple mapping patterns
This step has three parts: 1. realign with Bowtie2; 2. Filter realigned BAM; 3. Make REF and VAR depth table.

A warpper script is provided, which takes three positional arguments:
1. one miniBAM from **STEP 0**.
2. the identifier of the donor (e.g., patient1).
3. the output prefix (e.g., patient1_normal).
```bash
$ ./src/analyze_multi_map.sh $miniBAM $donor_id $out_prefix
```
Notes:
1. **Please modify the path to bowtie2 index in the warpper script**.
2. **Please make sure that the same patient identifier is used for each normal-tumour pair**.
3. The bash script also depends on other reference files as well as the location of the `src` directory if not in PATH.
4. No sanity checks are performed for inputs and references..

This step will generate multiple files:
1. `$out_prefix.realn.bam`: realigned BAM without filering
2. `$out_prefix.filtered.realn.bam`: filtered version of 1
3. `$out_prefix.unweighted.dp.tsv`: REF and VAR depth without multiple mapping adjustment
4. `$out_prefix.weighted.dp.tsv`: REF and VAR depth with multiple mapping adjustment
5. `$out_prefix.mstat.tsv`: multiple mapping pattern table

An example of each output file is provided in [./example/](./example/).

## STEP 2: Call somatic mutations
Our mutation calls depends on panel of normals (PON).
To make PON, we can simply combine all normal depth files from **STEP 1**.
For example, if all normal depth files are stored in the same directory with prefix like `patientID_normal`, we can do
```bash
# Get header
head -1 ./patient1_normal.weighted.dp.tsv > ../PON.weighted.dp.tsv
# Iterate through files
for f in ./*normal.weighted.dp.tsv
do
  tail -n +2 $f >> ../PON.weighted.dp.tsv
done
```
We can now call mutations with PON and tumour depth file (here using patient1 as an example):
```bash
$ ./src/beta_binomial_model.py ./patient1_tumour.weighted.dp.tsv ../PON.weighted.dp.tsv ./ref/ref_U1U11.map patient1_genotype.tsv
```
Notes:
1. The paired normal must present in the PON, because we are calling somatic mutations (in tumour but not in paired normal).
2. The [./ref/ref_U1U11.map](./ref/ref_U1U11.map) is a file mapping each nucleotide to a U1 genes (GRCh37).
3. In the output file ([patient1_genotype.tsv](./example/patient1_genotype.tsv)), the last column (`GT`) is the genotype.
4. The cut-off of EB scores (`EB_tumour` and `EB_normal` in the output file)
we used to call mutations uses PON with 30<n<50 (`num_p` and `num_n` in output file). If your sample size is different, you might need to find optimal cut-offs by yourself.  

## References
1. **Realignment using Bowtie2**: Langmead, B. & Salzberg, S. L. Fast gapped-read alignment with Bowtie 2. Nat. Methods 9, 357–359 (2012).
2. **Mutation call using EBCall**: Shiraishi, Y. et al. An empirical Bayesian framework for somatic mutation detection from cancer genome sequencing data. Nucleic Acids Res. 41, e89 (2013).

## Bug report
These scripts are not fully tested, so please start a new issue if you need help or find a bug.
