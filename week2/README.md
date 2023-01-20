# Variant Calling

Perform germline variant calling on a normal sample using `GATK HaplotypeCaller`.

#### Tutorial Data

Navigate to your clone of the MA5112 directory from last week, making sure you are on your own branch:

```console
git branch --all
```

```console
* your_branch
  main
  remotes/origin/HEAD -> origin/main
  remotes/origin/main
```

> Your branch has still not been published to the repository i.e the remote origin.

Synchronise the week2 tutorial to your branch by running:

```console
git pull origin main
```

The data files are available under `week2/data`.

***

Try pushing your branch to the remote origin:

```console
git add .
git commit 'first commit'
git push
```

#### Tutorial packages

I have created a container for the tutorial which you can access via Docker Hub:

```console
docker pull barryd237/week2
```

Packages of note in the container are:

```console
gatk4-4.3
samtools-1.11
bwa-0.7.17
```

#### Mounting tutorial data

To bind the tutorial data to the container, we will use the `-v` flag again.

Make sure you are in the MA5112 directory:

```console
docker run -it -v $(pwd):/files/ barryd237/week2
```

The MA5112 directory is now bound to the container, in the directory `/files/`. Simply `cd` into this directory to access the tutorial files:

```console
cd files/week2/data/
ls -la
```

***

As an alternative, you can create a conda envrionment using the following `environment.yml` file:

```yaml
name: week2
channels:
 - bioconda
dependencies:
 - bwa
 - gatk4
 - samtools
```

This bypasses the need to mount containers.

***

![germs](../docs/images/germline.png)

## Genome Index

> Such programs ... use a computational strategy known as ‘indexing’ to speed up their mapping algorithms. Like the index at the end of a book, an index of a large DNA sequence allows one to rapidly find shorter sequences embedded within it. [_How to map billions of short reads onto genomes_](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2836519/)

```console
bwa index genome.fasta
```

Indexing with `samtools` and `gatk` is less computationally intensive, it simply outputs the contigs & contig lengths present in the reference.

```console
samtools faidx genome.fasta
gatk CreateSequenceDictionary -R genome.fasta -O genome.dict
```

## Genomic Intervals

Instead of targeted exome panels, I will show you how to parellize computationally intensive steps via chromosome arms. The exact same principles apply if you were to use an exome capture kit (same flags, same file format (BED)).

```console
head -n 1 multi_intervals.bed > tmp.bed && mv tmp.bed chr21_2-23354000.bed
tail -n 1 multi_intervals.bed > tmp.bed && mv tmp.bed chr21_25689498-46709983.bed
```

Downstream variant calling tools demand files are compressed using `bgzip` & indexed using `tabix`:

```console
bgzip chr21_2-23354000.bed
tabix -p bed chr21_2-23354000.bed.gz
```

```console
bgzip chr21_25689498-46709983.bed
tabix -p bed chr21_25689498-46709983.bed.gz
```

## Genome Alignment

In the field of variant calling, `bwa` is the academic standard for aligning DNA reads to the reference genome.

BWA requires the indexed reference genome file, a read group identifier and the input fastq files. We will use a unix pipe `|` to pass the output of BWA directly to the `samtools sort` command to sort the genomic alignments by coordinates/names.

```console
bwa mem \
    genome.fasta \
    -R "@RG\tID:HCC1395\tSM:HCC1395N\tPL:ILLUMINA" \
    HCC1395N_R1.fastq.gz HCC1395N_R2.fastq.gz | samtools sort - > HCC1395N.bam
```

## Mark Duplicates

> “Almost all statistical models for variant calling assume some sort of independence between measurements. The duplicates (if one assumes that they arise from PCR artifact) are not independent. This lack of independence will usually lead to a breakdown of the statistical model and measures of statistical significance that are incorrect” – Sean Davis.

```console
gatk \
    --java-options -Xmx1g \
    MarkDuplicates \
    -I HCC1395N.bam \
    -M HCC1395N.mkdup.metrics \
    -O HCC1395N.mkdup.bam \
    --CREATE_INDEX true
```

## Base quality score recalibration

Base Quality Score Recalibration. A data pre-processing step that detects systematic errors made by the sequencing machine when it estimates the accuracy of each base call. IT also incorporates information about loci that are known to vary in populations using variant resource files (dbSNP, Mills KG).

:exclamation: We will split the analysis over the two chromosome arms of chr21

```console
gatk --java-options -Xmx1g \
     BaseRecalibrator \
     -I HCC1395N.mkdup.bam \
     -O HCC1395N.chr21_2-23354000.recal.table \
     -L chr21_2-23354000.bed.gz \
     -R genome.fasta \
     --known-sites mills_and_1000G.indels.hg38.vcf.gz \
     --known-sites dbsnp_138.hg38.vcf.gz
```

Run the same analysis again, using the second arm of chr21. You will need to adjust the `-O` and `-L` flags accordingly.

You should have two files when you are finished: `HCC1395N.chr21_2-23354000.recal.table` & `HCC1395N.chr21_25689498-46709983.recal.table`.

## Apply BQSR 

Apply the tables generated in the previous step to the BAM files to produce recalibrated (scores) BAM files.

```console
gatk --java-options -Xmx1g \
     ApplyBQSR \
     -I HCC1395N.mkdup.bam \
     -O HCC1395N.chr21_2-23354000.recal.bam \
     -L chr21_2-23354000.bed.gz \
     -R genome.fasta \
     --bqsr-recal-file HCC1395N.chr21_2-23354000.recal.table
```

:exclamation: Apply the above to the other arm in chromosome 21.

## Haplotype Caller

Call germline variants using `HaplotypeCaller`

```console
gatk --java-options -Xmx2g \
    HaplotypeCaller \
    -I HCC1395N.chr21_2-23354000.recal.bam \
    -O HCC1395N.haplotypecaller.chr21_2-23354000.vcf.gz \
    -R genome.fasta \
    -D dbsnp_138.hg38.vcf.gz \
    -L chr21_2-23354000.bed.gz
```

:exclamation: as before, carry out the same for the second arm of chr21.

## Combine VCFs

Combine the VCF files before carrying out variant filtering:

> If you used different output names for your second arm of chr21, apply it below.
> You need the two haplotypecaller.vcf files for this step...

```console
gatk --java-options -Xmx1g \
    MergeVcfs \
    -I HCC1395N.haplotypecaller.chr21_25689498-46709983.vcf.gz \
    -I HCC1395N.haplotypecaller.chr21_2-23354000.vcf.gz \
    -O HCC1395N.haplotypecaller.vcf.gz \
    --SEQUENCE_DICTIONARY genome.dict
```

## CNN score variants

Annotate a VCF with scores from a Convolutional Neural Network (CNN). This tool streams variants and their reference context to a python program, which evaluates a pre-trained neural network on each variant. The default models were trained on **single-sample VCFs**. The default model should not be used on VCFs with annotations from joint call-sets.

```console
gatk --java-options -Xmx1g \
    CNNScoreVariants \
    -V HCC1395N.haplotypecaller.vcf.gz \
    -O HCC1395N.cnn.vcf.gz \
    -R genome.fasta \
    -L multi_intervals.bed
```

> You will get an error about Python libraries, so I prepared one for you. The file is available under the `MA5112/week2/data/cnn/HCC1395N.cnn.vcf.gz` directory.

## Filter variants

Finally, filter the variants, incorporating information from known databases and CNN scoring.

```console
gatk --java-options -Xmx2g\
    FilterVariantTranches \
    -V cnn/HCC1395N.cnn.vcf.gz \
    --resource dbsnp_138.hg38.vcf.gz \
    --resource mills_and_1000G.indels.hg38.vcf.gz \
    --output HCC1395N.haplotypecaller.filtered.vcf.gz \
    --info-key CNN_1D \
    --indel-tranche 0
```