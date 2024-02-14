#!/bin/bash

# Function to remove file extension
remove_extension() {
    # Strip the extension from the input file to use as the default output prefix
    local filename=$(basename -- "$1")
    local name="${filename%.*}"
    echo "$name"
}

# Check for number of arguments
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <input_23andMe_file> [<output_prefix>]"
    exit 1
fi

# Input file is always the first argument
INPUT_FILE=$1

# If an output prefix isn't provided, use the input filename
if [ "$#" -eq 1 ]; then
    OUTPUT_PREFIX=$(remove_extension "$INPUT_FILE")
else
    OUTPUT_PREFIX=$2
fi

PYTHON_SCRIPT="23andme_to_ped.py" # Your provided python script

# Step 1: Convert 23andMe format to PED and MAP using 23andMe_to_PLINK.py
python3 $PYTHON_SCRIPT $INPUT_FILE $OUTPUT_PREFIX

# Step 2: Use PLINK2 to convert PED/MAP files to BED/BIM/FAM format
plink --file $OUTPUT_PREFIX --make-bed --out $OUTPUT_PREFIX

# Check if the conversion was successful
if [ $? -eq 0 ]; then
    echo "Conversion to PLINK BED format completed successfully."
else
    echo "Conversion to PLINK BED format failed."
    exit 1
fi