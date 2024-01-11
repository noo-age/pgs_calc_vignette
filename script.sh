#!/bin/bash
#bin of generally useful functions

#check target build before using
nextflow run pgscatalog/pgsc_calc \
    -profile docker \
    --input samplesheet.csv \
    --target_build GRCh37 \
    --pgs_id PGS000758 \


#recodes ped to bed
plink --file <ped_input_files> --make-bed --out <bed_output>

#scores using plink --score
plink2 --bfile <bed_input_files> --score <score_file> 1 4 6 --out <outfile>
plink2 --pfile GRCh37_1000G_ALL --score score_files/ea.txt 1 4 6 --out 1kgenomes