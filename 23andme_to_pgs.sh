#!/bin/bash

# Function to remove file extension
remove_extension() {
    # Strip the extension from the input file to use as the default output prefix
    local filename=$(basename -- "$1")
    local name="${filename%.*}"
    echo "$name"
}

# Check for the number of arguments
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input_23andMe_file> <score_file> [<output_prefix>]"
    exit 1
fi

# Input file and score file are always the first two arguments
INPUT_FILE=$1
SCORE_FILE=$2

# If an output prefix isn't provided, use the input filename
if [ "$#" -eq 2 ]; then
    OUTPUT_PREFIX=$(remove_extension "$INPUT_FILE")
else
    OUTPUT_PREFIX=$3
fi

PYTHON_SCRIPT="23andme_to_ped.py" # Your provided python script

# Step 1: Convert 23andMe format to PED and MAP using 23andMe_to_PLINK.py
python3 $PYTHON_SCRIPT $INPUT_FILE $OUTPUT_PREFIX

# Step 2: Use PLINK2 to convert PED/MAP files to BED/BIM/FAM format
plink --file $OUTPUT_PREFIX --make-bed --out $OUTPUT_PREFIX

# Step 3: Use PLINK2 to provide a Polygenic Score (PGS) of the input genotype
plink2 --bfile $OUTPUT_PREFIX --score $SCORE_FILE 1 4 6 no-mean-imputation --out ${OUTPUT_PREFIX}_pgs

# Check if the PGS calculation was successful
if [ $? -eq 0 ]; then
    echo "PGS calculation completed successfully, output can be found in ${OUTPUT_PREFIX}_pgs.sscore"
else
    echo "PGS calculation failed."
    exit 1
fi