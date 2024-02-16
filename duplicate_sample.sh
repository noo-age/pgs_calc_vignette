#!/bin/bash

# Check if correct number of arguments provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input.vcf.gz>"
    exit 1
fi

# Input and output from command line arguments
input_vcf=$(basename $1)
echo $input_vcf
input_dir=$(dirname $1)/
echo $input_dir

# Check to see if duplication is necessary
num_samples=$(bcftools query -l "${input_dir}${input_vcf}" | wc -l)
if [ "$num_samples" -gt 20 ]; then
    echo "Input VCF has more than 20 samples. Duplication unnecessary to use the Michigan imputation server."
    #exit 0
fi

# Get the sample name from the input VCF
sample=$(bcftools query -l "${input_dir}${input_vcf}" | head -n 1)

# An array to store the names of the temporary VCF files
temp_files=()

# Create 20 copies of the VCF with different sample names
for i in $(seq 1 30); do
    # Create a new sample name for this copy
    new_sample="${sample}_copy_${i}"
    
    # Create a new temporary VCF file name for this copy
    temp_file="temp_${i}.vcf.gz"
    temp_files+=("$temp_file")
    
    # Use bcftools reheader to change the sample name
    bcftools reheader -s <(echo "$new_sample") -o "$temp_file" "${input_dir}${input_vcf}"
    bcftools index "$temp_file"
done

# Merge the temporary VCF files into the output VCF
bcftools merge "${temp_files[@]}" -f -Oz -o "${input_dir}dup-${input_vcf}"
bcftools index ${input_dir}dup-${input_vcf}

# Clean up the temporary files
for temp_file in "${temp_files[@]}"; do
    rm "$temp_file"
    rm "${temp_file}.csi" 
done

echo "Sample duplication complete"