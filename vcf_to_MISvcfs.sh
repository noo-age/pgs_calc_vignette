#!/bin/bash

# Check if an input file and outdir were provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.vcf> <outdir>"
    exit 1
fi

input_vcf=$1
outdir=$2

# bgzip input if it isn't already
extension="${input_vcf##*.}"
if [ "$extension" != "gz" ]; then
    bgzip -c $input_vcf > ${input_vcf}.gz
    input_vcf=${input_vcf}.gz
fi

input_vcf=$(basename $1)
input_dir=$(dirname $1)/

chromosome_split_script="./filter_split_vcf.sh" 
duplicate_sample_script="./duplicate_sample.sh"

# Create a duplicate of the input VCF
$duplicate_sample_script ${input_dir}/${input_vcf}

# Split the duplicated VCF by chromosomes and output to the specified outdir
$chromosome_split_script ${input_dir}/dup-${input_vcf} $outdir

