# ChIP-Seq

Training materials generously provided by galaxy project.. 

![chip](../docs/images/chip-seq.png)

## Analysis BAM files

In the interest of time, BAM files have been prepared for you by aligning paired-end DNA sequencing reads to the reference genome.

* Create a directory `bams/` to store the BAM files.

* Each entry in the `files.txt` file contains the URL for an analysis BAM file. Loop over the contents of the file, applying the `wget` command to download them to your working directory. (_hint: use a while loop_).

***

Recall last week we used `samtools` to view the header of the BAM files.

1. What tool was used to align reads to the reference genome?

2. What processing was performed on the BAM file?

***

* Merge the BAM replicates using `samtools merge`.

* Index the files (`samtools index`) and delete the original replicate BAM files.


## Quality Control

Recall that ChIP-Seq experiments require a background control sample - the 'input' sample - one that has been cross-linked and sonicated **but not immuno-precipitated**.

We expect our input sample to have a flat distribution across the sequenced genome whilst immuno-precipitated proteins of interest will have pileups - 'peaks' - of reads in the regions in which they were bound to DNA.

* Make a directory `deeptools/qc` to store the output files from QC.

### Sample Heterogeneity

We will assess sample heterogeneity by computing the correlation of read counts on different regions for all samples.

> This is analogous to RNA-Seq which we will cover next week. Instead of counting reads that span a transcript, ChIP-Seq involves counting reads that span a region.

We expect that the replicates of the ChIP-seq experiments should be clustered more closely to each other than the replicates of the input sample.

### multiBamSummary

Compare the number of reads spanning a pre-defined bin range (`--binSize`) between samples. We can reduce the computational load by sampling reads in 1Kb bins every 500bp intervals (`--distanceBetweenBins`).

```console
multiBamSummary bins \
    --outFileName deeptools/multisummary.matrix \
    --binSize 1000 \
    --distanceBetweenBins 500 \
    --region chrX \
    --bamfiles bams/wt_CTCF.bam bams/wt_H3K27me3.bam bams/wt_H3K4me3.bam bams/wt_input.bam
```

The output file is compressed, we cannot view it using unix tools. 

### PlotCorrelation

Plot a sample-sample heatmap using the `plotCorrelation` tool:

```console
plotCorrelation \
    -in deeptools/multisummary.matrix \
    --whatToPlot heatmap \
    --corMethod pearson \
    -o deeptools/multisummary_pearson.pdf
```

* Take note of the output from the terminal and apply the suggested changes to a new pdf file.

***

1. Why are the suggested changes important when using `--corMethod pearson` ?
2. 


#### plotPCA 

```console
plotPCA \
    --corData deeptools/multisummary.matrix \
    --transpose \
    --plotFile deeptools/multisample_pca.pdf
```

#### plotFingerprint


```console
plotFingerprint \
    -b bams/wt_CTCF.bam bams/wt_H3K27me3.bam bams/wt_H3K4me3.bam bams/wt_input.bam \
    -plot deeptools/fingerprint.pdf \
    --region chrX \
    --numberOfSamples 10000
```

failed samples or broad peaks? domain knowledge is important here. slighly worrying its considered broader than input though. 

## BedGraph

### normalisation

#### bedgraph

bamCoverage --bam bams/wt_input.bam --binSize 25 --region chrX --outFileName deeptools/input.bedgraph --outFileFormat bedgraph --effectiveGenomeSize 2308125349 --normalizeUsing RPGC


bamCoverage --bam bams/wt_H3K4me3.bam --binSize 25 --region chrX --outFileName deeptools/h3k4me.bedgraph --outFileFormat bedgraph --effectiveGenomeSize 2308125349 --normalizeUsing RPGC

which regions in the bedgraph file have the highest coverage? hint: use unix sort on the fourth column, in descending order. google is your friend here. 

run the same for the BigWig - goal is to get you familiar with these file formats. 

#### bigwig

bamCoverage --bam bams/wt_input.bam --binSize 25 --region chrX --outFileName deeptools/input.bigwig --outFileFormat bigwig --effectiveGenomeSize 2308125349 --normalizeUsing RPGC


bamCoverage --bam bams/wt_H3K4me3.bam --binSize 25 --region chrX --outFileName deeptools/h3k4me.bigwig --outFileFormat bigwig --effectiveGenomeSize 2308125349 --normalizeUsing RPGC


## normalise h3k4me using bamcompare

normalise against the input control to isolate true ChIP-Seq peaks 
 ORDER matters here! --bamfile1 is the treated sample

```console
bamCompare --bamfile1 bams/wt_H3K4me3.bam --bamfile2 bams/wt_input.bam --region chrX --outFileName deeptools/H3K4me_vs_input.bigwig --outFileFormat bigwig --operation log2
```

autoscale your IGV browser or else its muck. poay attention to the scale [0, number] in IGV. 


macs yo

make macs2 directory

macs2 callpeak -t bams/wt_H3K4me3.bam -c bams/wt_input.bam -f BAMPE -g mm -n H3K4me3_experiment --outdir macs2/


## plotting vs two samples

the steps computed above should be done for all samples vs. input. do input vs CTCF now

bamCompare --bamfile1 bams/wt_CTCF.bam --bamfile2 bams/wt_input.bam --region chrX --outFileName deeptools/CTCF_vs_input.bigwig --outFileFormat bigwig --operation log2

macs2 callpeak -t bams/wt_CTCF.bam -c bams/wt_input.bam -f BAMPE -g mm -n CTCF_experiment --outdir macs2/

cat macs2/H3K4me_experiment_summits.bed macs2/CTCF_experiment_summits.bed | awk -vFS="\t" '{print $1,$2,$3}' - | sort -k1,1 -k2,2n | uniq | tr ' ' '\t'  > H3K4me_CTCF_merged.bed



use two vs bigwig files from bamcompare

computeMatrix reference-point --regionsFileName H3K4me3_CTCF_merged.bed --scoreFileName deeptools/H3K4me3_vs_input.bigwig deeptools/CTCF_vs_input.bigwig --upstream 3000 --downstream 3000 --referencePoint center --outFileName H3K4me3_CTCF.matrix
