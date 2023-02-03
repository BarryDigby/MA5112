#!/bin/bash

name: week1
channels:
  - bioconda
dependencies:
  - fastqc
  - multiqc
  - trim-galore


read -p "Enter path for data for adapter/read trimming:" p1
read -p "Enter desired output directory name:" p2

mkdir -p p1/p2 +"/raw"

for file in $p1/*; do fastqc $file --outdir p2; done


mkdir -p 
