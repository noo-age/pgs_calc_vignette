#!/bin/bash

# Function to remove file extension
remove_extension() {
    local filename=$(basename -- "$1")
    local name="${filename%.*}"
    echo "$name"
}

# Check for the correct number of arguments
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <filetype> <input_file> <score_file> [<output_prefix>]"
    exit 1
fi

# Filetype and input file are the first two arguments
FILETYPE=$1
INPUT_FILE=$2
SCORE_FILE=$3

# If an output prefix isn't provided, use the input filename
if [ "$#" -eq 3 ]; then
    OUTPUT_PREFIX=$(remove_extension "$INPUT_FILE")
else
    OUTPUT_PREFIX=$4
fi

# Processing based on filetype
case $FILETYPE in
    23andme)
        PYTHON_SCRIPT="23andme_to_ped.py"
        python3 $PYTHON_SCRIPT $INPUT_FILE $OUTPUT_PREFIX
        plink --file $OUTPUT_PREFIX --make-bed --out $OUTPUT_PREFIX
        ;;
    bed)
        # Assuming BED/BIM/FAM files exist with $OUTPUT_PREFIX, proceed directly to scoring
        ;;
    vcf)
        plink2 --vcf $INPUT_FILE --make-bed --out $OUTPUT_PREFIX
        ;;
    *)
        echo "Invalid filetype specified. Valid options are: 23andme, bed, vcf"
        exit 1
        ;;
esac

# Perform Polygenic Score (PGS) analysis
plink2 --bfile $OUTPUT_PREFIX --score $SCORE_FILE 1 4 6 no-mean-imputation --out ${OUTPUT_PREFIX}_pgs

# Check if the PGS calculation was successful
if [ $? -eq 0 ]; then
    echo "PGS calculation completed successfully, results are in: ${OUTPUT_PREFIX}_pgs.sscore"
else
    echo "PGS calculation failed."
    exit 1
fi