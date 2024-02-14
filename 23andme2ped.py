import pandas as pd
import sys

# This function parses the 23andMe raw data file and extracts the necessary fields
def parse_23andMe(file_path):
    with open(file_path, 'r') as file:
        lines = [line.strip() for line in file if not line.startswith('#') and line.strip()] # Remove comments and empty lines
    
    # Split each line by tabs and create tuples (rsid, chromosome, position, genotype)
    records = [line.split('\t') for line in lines if line] # Remove any blank lines that are left
    
    # Filter out X and Y chromosome data
    filtered_records = [(record[0], record[1], record[2], record[3]) for record in records if record[1].upper() not in ['X', 'Y', 'MT']]
    
    if not filtered_records:
        return [], [], [], []

    rsids, chromosomes, positions, genotypes = zip(*filtered_records)
    
    return rsids, chromosomes, positions, genotypes

# This function writes the .ped and .map files which are necessary for PLINK
def write_ped_map(rsids, chromosomes, positions, genotypes, output_prefix):
    # Create .map file content
    map_content = [f'{chrom}\t{rsid}\t0\t{pos}' for rsid, chrom, pos in zip(rsids, chromosomes, positions)]

    # Split genotypes from 'AA' to 'A A' for PLINK .ped format to read
    split_genotypes = [' '.join(list(genotype)) for genotype in genotypes]

    # Create .ped file content. We assume the individual's ID, family ID, etc., as PLINK requires it.
    # Here we use dummy values -- you should replace 'FID', 'IID', 'PID', 'MID', 'Sex' and 'Phenotype' with actual information if known.
    ped_content = ["FID IID PID MID Sex Phenotype " + ' '.join(split_genotypes)]

    # Write .map file
    with open(f'{output_prefix}.map', 'w') as map_file:
        map_file.write('\n'.join(map_content))

    # Write .ped file
    with open(f'{output_prefix}.ped', 'w') as ped_file:
        ped_file.write('\n'.join(ped_content))


def main():
    # Check for command line arguments for the input file and output prefix
    if len(sys.argv) < 3:
        print('Usage: python3 23andMe_to_PLINK.py <input_23andMe_file> <output_prefix>')
        sys.exit(1)

    input_file = sys.argv[1]
    output_prefix = sys.argv[2]
    
    # Parse the 23andMe file
    rsids, chromosomes, positions, genotypes = parse_23andMe(input_file)
    # Convert to .ped and .map format
    write_ped_map(rsids, chromosomes, positions, genotypes, output_prefix)

if __name__ == '__main__':
    main()