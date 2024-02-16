#!/bin/bash

# Check if an input file and output directory were provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_vcf.gz> <outdir>"
    exit 1
fi

# The input VCF file from the first command line argument
input_vcf=$(basename $1)
input_dir=$(dirname $1)

# The output directory from the second command line argument
outdir="$2"

# Create the output directory if it doesn't exist
mkdir -p "$outdir"

# Loop through chromosomes 1-22
for i in {1..22} X Y
do
    # Define chromosome name
    chr="chr${i}"

    # Define output file name
    output_vcf="${outdir}/${chr}_${input_vcf}"
    
    # Filter the input VCF for the current chromosome and compress the output
    bcftools view ${input_dir}/${input_vcf} --regions "$chr" | bgzip > "$output_vcf"
    
    # Index the output VCF file
    bcftools index "$output_vcf"
done

echo "Chromosome VCF files have been generated in $outdir."