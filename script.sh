nextflow run pgscatalog/pgsc_calc \
    -profile docker \
    --input samplesheet.csv \
    --target_build GRCh37 \
    --pgs_id PGS001229 \


#recodes ped to bed
plink --file <ped_input_files> --make-bed --out <bed_output>

#scores using plink --score
plink2 --bfile <bed_input_files> --score <score_file> 1 4 6 --out <outfile>