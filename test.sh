#!/usr/bin/env bash

## Kallisto indexing step

data_dir="/data/github/MA5112/week4/toy_data"

kallisto index -i ${data_dir}/test.idx ${data_dir}/transcripts.fasta

## Kallisto quant step

kallisto quant -i ${data_dir}/test.idx -o ${data_dir}/results ${data_dir}/reads_1.fastq.gz ${data_dir}/reads_2.fastq.gz
