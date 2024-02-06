#!/bin/bash

# Check if correct number of arguments provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_vcf> <output_vcf>"
    exit 1
fi

# Input and output from command line arguments
input_vcf="$1"
output_vcf="$2"

# Get the sample name from the input VCF
sample=$(bcftools query -l "$input_vcf")

# An array to store the names of the temporary VCF files
temp_files=()

# Create 20 copies of the VCF with different sample names
for i in $(seq 1 20); do
    # Create a new sample name for this copy
    new_sample="${sample}_copy_${i}"
    
    # Create a new temporary VCF file name for this copy
    temp_file="temp_${i}.vcf.gz"
    temp_files+=("$temp_file")
    
    # Use bcftools reheader to change the sample name
    bcftools reheader -s <(echo "$new_sample") -o "$temp_file" "$input_vcf"
    bcftools index "$temp_file"
done

# Merge the temporary VCF files into the output VCF
bcftools merge "${temp_files[@]}" -Oz -o "$output_vcf"
bcftools index "$output_vcf"

# Clean up the temporary VCF files
rm "${temp_files[@]}"

echo "Finished duplicating the sample in the VCF file."