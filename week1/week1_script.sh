#!/usr/bin/env bash 

# set path to analysis files
DATA=$1

# 1. Perform QC on raw data

mkdir -p quality_control
mkdir -p quality_control/raw

for file in ${DATA}/*; do
    fastqc $file \
        --outdir quality_control/raw
done

multiqc quality_control/raw -n raw_reads -o quality_control

# 2. Trim and QC

mkdir -p quality_control/trimmed

for file in ${DATA}/*; do

    trim_galore \
        --adapter TGGAATTCTCGGGTGCCAAGG \
        --length 17 \
        --clip_r1 4 \
        --three_prime_clip_r1 4 \
        --max_length 30 \
        --gzip \
        --fastqc \
        --fastqc_args "--outdir quality_control/trimmed" \
        --output_dir quality_control/trimmed \
        $file

done

multiqc quality_control/trimmed -n trimmed_reads -o quality_control
