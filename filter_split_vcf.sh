#!/bin/bash

# Check if an input file was provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_vcf.gz>"
    exit 1
fi

# The input VCF file from the first command line argument
input_vcf="$1"

# Loop through chromosomes 1-22
for i in {1..22}
do
    # Define chromosome name
    chr="chr${i}"
    
    # Extract the base name of the input file without its extension
    base_name=$(basename "$input_vcf" .vcf.gz)
    
    # Define output file name
    output_vcf="${chr}_${base_name}.vcf.gz"
    
    # Filter the input VCF for the current chromosome and compress the output
    bcftools view "$input_vcf" --regions "$chr" | bgzip > "$output_vcf"
    
    # Index the output VCF file
    tabix -p vcf "$output_vcf"
done

echo "Chromosome VCF files have been generated."
