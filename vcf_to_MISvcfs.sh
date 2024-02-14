#!/bin/bash

# Check if an input file and outdir were provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.vcf.gz> <outdir>"
    exit 1
fi

input_vcf="$1"
outdir="$2"
chromosome_split_script="./filter_split_vcf.sh" 
duplicate_sample_script="./duplicate_sample.sh"

# Create a duplicate of the input VCF
$duplicate_sample_script $input_vcf $input_vcf

# Split the duplicated VCF by chromosomes and output to the specified outdir
$chromosome_split_script $input_vcf $outdir

