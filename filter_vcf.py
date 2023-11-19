import sys

def is_standard_chromosome(chrom):
    # Define standard chromosome names (1-22, X, Y, MT)
    standard_chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y', 'MT']
    return chrom.lstrip('chr') in standard_chromosomes

def filter_vcf(input_vcf_path):
    # Generate the output VCF path by adding 'filtered' before the file's extension
    output_vcf_path = input_vcf_path.rsplit('.', 1)[0] + '_filtered.vcf'
    with open(input_vcf_path, 'r') as input_vcf, open(output_vcf_path, 'w') as output_vcf:
        for line in input_vcf:
            line = line.strip()
            # Write header lines directly to the output file
            if line.startswith('#'):
                output_vcf.write(line + '\n')
            else:
                # Process variant lines
                parts = line.split('\t')
                chrom = parts[0]
                alt_alleles = parts[4].split(',')
                # Filter for standard chromosomes and biallelic variants
                if is_standard_chromosome(chrom) and len(alt_alleles) == 1:
                    output_vcf.write(line + '\n')
    print(f"Filtered VCF file created: {output_vcf_path}")

# Check if a command-line argument was provided for the input VCF file path
if len(sys.argv) != 2:
    print("Usage: python3 filter_vcf.py <input_vcf_path>")
    sys.exit(1)

input_vcf_file_path = sys.argv[1]

filter_vcf(input_vcf_file_path)