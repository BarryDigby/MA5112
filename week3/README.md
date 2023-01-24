# ChIP-Seq

Training materials generously provided by galaxy project...

Recommend using the desktop computers in ADB for their large screen size when using IGV.

## Analysis BAM files

In the interest of time, BAM files have been prepared for you by aligning paired-end DNA sequencing reads to the reference genome.

* Create a directory `bams/` to store the BAM files.

* Each entry in the `files.txt` file contains the URL for an analysis BAM file. Loop over the contents of the file, applying the `wget` command to download them to your working directory. (_hint: use a while loop_).

* Index the files (`samtools index`).

***

Recall last week we used `samtools` to view the header of the BAM files.

1. What tool was used to align reads to the reference genome?

2. What processing was performed on the BAM file?

## Quality Control

Recall that ChIP-Seq experiments require a background control sample - the 'input' sample - one that has been cross-linked and sonicated **but not immuno-precipitated**.

We expect our input sample to have a flat distribution across the sequenced genome whilst immuno-precipitated proteins of interest will have pileups - 'peaks' - of reads in the regions in which they were bound to DNA.

* Make a directory `deeptools/` to store the output files from QC.

### Sample Heterogeneity

We will assess sample heterogeneity by computing the correlation of read counts on different regions for all samples.

> This is analogous to RNA-Seq which we will cover next week. Instead of counting reads that span a transcript, ChIP-Seq involves counting reads that span a region.

We expect that the replicates of the ChIP-seq experiments should be clustered more closely to each other.

* Compare the number of reads spanning a pre-defined bin range (`--binSize`) between samples using `multiBamSummary bins`. We can reduce the computational load by sampling reads in 1Kb bins every 500bp intervals (`--distanceBetweenBins`).

```console
multiBamSummary bins \
    --outFileName deeptools/multisummary.matrix \
    --binSize 1000 \
    --distanceBetweenBins 500 \
    --region chrX \
    --bamfiles bams/wt_CTCF_rep1.bam bams/wt_CTCF_rep2.bam bams/wt_H3K27me3_rep1.bam bams/wt_H3K27me3_rep2.bam bams/wt_H3K4me3_rep1.bam bams/wt_H3K4me3_rep2.bam bams/wt_input_rep1.bam bams/wt_input_rep2.bam
```

* Plot a sample-sample heatmap using the `plotCorrelation` tool:

```console
plotCorrelation \
    -in deeptools/multisummary.matrix \
    --whatToPlot heatmap \
    --corMethod pearson \
    -o deeptools/multisummary_pearson.pdf
```

* Take note of the output from the terminal and apply the suggested changes to a new pdf file.

1. Why are the suggested changes important when using `--corMethod pearson` ?
2. Is the heatmap indicative of high quality experimental design?

***

Create a principal component analysis plot between PC1 vs. PC2:

```console
plotPCA \
    --corData deeptools/multisummary.matrix \
    --transpose \
    --plotFile deeptools/multisample_pca.pdf
```

1. Do replicates cluster closely together?
2. Which ChIP-Seq sample is closely related to the input control? In terms of peaks, what do you think this might mean?

### IP Strength

An ideal ‘input’ with perfect uniform distribution of reads along the genome (i.e. without enrichments in open chromatin) and infinite sequencing coverage should generate a straight diagonal line. A very specific and strong ChIP enrichment will be indicated by a prominent and steep rise of the cumulative sum towards the highest rank. This means that a big chunk of reads from the ChIP sample is located in few bins which corresponds to high, narrow enrichments typically seen for transcription factors.

If your ChIP-Seq sample is close to the input line, this does not indicate a failed experiment - just broader peaks!

```console
plotFingerprint \
    -b bams/wt_H3K4me3_rep1.bam bams/wt_input_rep1.bam \
    -plot deeptools/fingerprint.pdf \
    --region chrX \
    --numberOfSamples 10000
```

Read how to interpret the fingerprint plot here: [deeptools: fingerprint plot](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html#what-the-plots-tell-you).

***

* Now that we are satisfied our replicates are of decent quality, merge the BAM replicates using `samtools merge`.

* Index the merged bam files.

## Sample Normalisation

### Motivation

Generate BedGraph coverage files of the `H3Kme3` and `Input` samples:

```console
bamCoverage --bam bams/wt_input.bam --binSize 25 --region chrX --outFileName deeptools/input.bedgraph --outFileFormat bedgraph --effectiveGenomeSize 2308125349 --normalizeUsing RPGC
bamCoverage --bam bams/wt_H3K4me3.bam --binSize 25 --region chrX --outFileName deeptools/h3k4me3.bedgraph --outFileFormat bedgraph --effectiveGenomeSize 2308125349 --normalizeUsing RPGC
```

1. What do each of the columns in the bedgraph file represent?
2. which regions in each bedgraph file have the highest coverage? (_hint: use unix sort on the fourth column, in descending order. google is your friend here..._)

Generate the same coverage using bigwig files

```console
bamCoverage --bam bams/wt_input.bam --binSize 25 --region chrX --outFileName deeptools/input.bigwig --outFileFormat bigwig --effectiveGenomeSize 2308125349 --normalizeUsing RPGC
bamCoverage --bam bams/wt_H3K4me3.bam --binSize 25 --region chrX --outFileName deeptools/h3k4me.bigwig --outFileFormat bigwig --effectiveGenomeSize 2308125349 --normalizeUsing RPGC
```

Some notes on BedGraph vs BigWig are available here: [graphing track data formats](http://genomewiki.ucsc.edu/index.php/Selecting_a_graphing_track_data_format)

#### IGV

Open a new terminal and run IGV (`conda activate week3 && igv`). Load the following files in IGV: `input.bigwig`, `wt_input.bam`, `h3k4me3.bigwig`, `h3k4me3.bam`.

* Remove the BAM coverage track.
* Navigate to `chrX:48,249,729-48,267,236`
* Right click the bigwig files and AutoScale the data

Hopefully you can see the even distribution of reads in the input sample vs the H3K4me3 sample. The next step involves normalising the H3K4me3 sample against the input for peak callin gusing `macs2`.

## Normalisation

Normalise the H3K4me3 sample using the input control:

```console
bamCompare --bamfile1 bams/wt_H3K4me3.bam --bamfile2 bams/wt_input.bam --region chrX --outFileName deeptools/H3K4me_vs_input.bigwig --outFileFormat bigwig --operation log2
```

Load the output file `H3K4me_vs_input.bigwig` into your running IGV session. Be sure to Autoscale the data.

## Detect Enriched Peaks

We could see areas of enrichment in IGV. Use `macs2` to detect and score these peaks computationally:

* Make a directory `mcas2/`

```console
macs2 callpeak -t bams/wt_H3K4me3.bam -c bams/wt_input.bam -f BAMPE -g mm -n H3K4me3_experiment --outdir macs2/
```

* Load both the `*_summits.bed` and `*_peaks.narrowPeak` file into IGV. What is the difference between the two files?

## Heatmaps vs two samples

Run through the steps below in your own time to plot a heatmap vs `H3K4me3` & `CTCF`:

* Normalise the `CTCF` sample using `input`:

```console
bamCompare --bamfile1 bams/wt_CTCF.bam --bamfile2 bams/wt_input.bam --region chrX --outFileName deeptools/CTCF_vs_input.bigwig --outFileFormat bigwig --operation log2
```

* Detect enriched peaks using macs2:


```console
macs2 callpeak -t bams/wt_CTCF.bam -c bams/wt_input.bam -f BAMPE -g mm -n CTCF_experiment --outdir macs2/
```

* Concatenate the outputs from macs2 and sort them, we only need Chr, start, end positions:

```console
cat macs2/H3K4me_experiment_summits.bed macs2/CTCF_experiment_summits.bed | awk -vFS="\t" '{print $1,$2,$3}' - | sort -k1,1 -k2,2n | uniq | tr ' ' '\t'  > H3K4me_CTCF_merged.bed
```

`ComputeMatrix` calculates scores per genome regions and prepares an intermediate file that can be used with plotHeatmap and plotProfiles. Typically, the genome regions are genes, but any other regions defined in a BED file can be used. `computeMatrix` accepts multiple score files (bigWig format) and multiple regions files (BED format). This tool can also be used to filter and sort regions according to their score.

```console
computeMatrix reference-point --regionsFileName H3K4me3_CTCF_merged.bed --scoreFileName deeptools/H3K4me3_vs_input.bigwig deeptools/CTCF_vs_input.bigwig --upstream 3000 --downstream 3000 --referencePoint center --outFileName H3K4me3_CTCF.matrix
```

* Generate a sample vs sample heatmap (all regions are compared)

```console
plotHeatmap --matrixFile H3K4me_CTCF.matrix --outFileName H3K4me_CTCF_2.pdf --kmeans 2
```