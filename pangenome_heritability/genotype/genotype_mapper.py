import os
import click
import pandas as pd
import re
import subprocess
import pysam
from typing import Dict, NamedTuple
from ..config import Config  # Ensure this module contains the Config class

class PlinkFiles(NamedTuple):
    bed: str
    bim: str
    fam: str

def load_csv(file_path: str) -> pd.DataFrame:
    """Load the CSV file with explicit handling for column names."""
    try:
        csv_data = pd.read_csv(file_path, sep=',', encoding='utf-8')
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        raise

    csv_data.columns = csv_data.columns.str.strip()
    return csv_data

def parse_fasta(file_path: str) -> Dict[str, list]:
    """Parse the FASTA file to extract variant names and their sequences."""
    variants = {}
    current_group = None
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:]
                if "Variant" in header:
                    variant_name = header
                    if current_group:
                        variants.setdefault(current_group, []).append(variant_name)
                else:
                    current_group = header
    return variants

def replace_seq_with_variants(csv_data: pd.DataFrame, variants: Dict[str, list]) -> pd.DataFrame:
    """Replace 'sequence_id' with actual variant names from the FASTA file."""
    for index, row in csv_data.iterrows():
        group = row['chromosome_group']
        seq_id_str = row['sequence_id']
        if 'seq' in seq_id_str:
            seq_id = int(seq_id_str.replace('seq', '')) - 1
            if group in variants and seq_id < len(variants[group]):
                csv_data.at[index, 'sequence_id'] = variants[group][seq_id]
            else:
                print(f"No variant available for {group} at sequence {seq_id}")
    return csv_data

def process_comparison_column(csv_data: pd.DataFrame) -> pd.DataFrame:
    """Process the 'comparison' column to extract genotype information."""
    print("Columns before processing 'comparison':", csv_data.columns.tolist())
    
    if 'comparison' not in csv_data.columns:
        print("Warning: 'comparison' column is missing or empty.")
        return csv_data

    def process_comparison(s):
        try:
            cleaned_s = re.sub(r',\s*$', '', s)
            return eval(cleaned_s)
        except Exception as e:
            print(f"Error processing comparison data: {s}. Error: {e}")
            return None

    csv_data['comparison'] = csv_data['comparison'].apply(process_comparison)
    csv_data = csv_data[csv_data['comparison'].notnull()]
    return csv_data

def extract_group_from_chromosome_group(csv_data: pd.DataFrame) -> pd.DataFrame:
    """Extract 'group' information from 'chromosome_group' column."""
    if 'group' not in csv_data.columns:
        print("Generating 'group' column from 'chromosome_group'.")
        # Split based on the last underscore and take the last part (e.g., '2' from 'Group_2_1')
        csv_data['group'] = csv_data['chromosome_group'].apply(lambda x: x.split('_')[-1])
    return csv_data

def add_start_column_if_missing(csv_data: pd.DataFrame) -> pd.DataFrame:
    """Generate 'start' column from 'sequence_id'."""
    if 'start' not in csv_data.columns:
        print("Generating 'start' column from 'sequence_id'.")
        def extract_start(seq_id):
            match = re.match(r'Variant_(\d+)_(\d+)_(\d+)_(\d+)', seq_id)
            if match:
                return int(match.group(3))  # Extract start position
            else:
                return None
        csv_data['start'] = csv_data['sequence_id'].apply(extract_start)
        if csv_data['start'].isnull().any():
            raise ValueError("Failed to extract 'start' from some 'sequence_id' entries.")
    return csv_data

def create_ped_and_map_files(csv_data: pd.DataFrame, vcf_file: str, output_prefix: str):
    """Generate PED and MAP files from CSV and VCF data."""
    vcf = pysam.VariantFile(vcf_file)
    samples = list(vcf.header.samples)
    
    # Extract variant data from VCF
    variant_data = {}
    for record in vcf:
        chromosome = str(record.contig)  # Convert chromosome to string for key consistency
        start = str(record.pos)          # Convert start to string for key consistency
        key = (chromosome, start)
        variant_data[key] = {sample: record.samples[sample]['GT'] for sample in samples}
    vcf.close()

    # Combine genotype function
    def combine_genotypes(existing, new):
        def combine_gt(gt1, gt2):
            if None in gt1 or -1 in gt1:
                return gt2
            if None in gt2 or -1 in gt2:
                return gt1
            return tuple(max(a, b) for a, b in zip(gt1, gt2))
        return {sample: combine_gt(existing.get(sample, (0, 0)), new.get(sample, (0, 0))) for sample in samples}

    # Initialize PED and MAP data
    ped_data = {sample: [sample, sample, '0', '0', '0', '-9'] for sample in samples}
    map_data = []

    # Group by chromosome and group
    groups = csv_data.groupby(['chromosome', 'group'])

    # Process each group
    for (chromosome, group), group_data in groups:
       
        comparison_matrix = group_data['comparison'].tolist()
        if not comparison_matrix:
            print(f"Warning: Empty comparison matrix for group {group}.")
            continue

        num_kmers = len(comparison_matrix[0]) if comparison_matrix else 0

        for idx in range(num_kmers):
            merged_genotypes = {sample: (None, None) for sample in samples}

            # Merge genotype data
            for row_idx, comparison in enumerate(comparison_matrix):
                value = comparison[idx]

                if value == 1:
                    key = (str(group_data.iloc[row_idx]['chromosome']), str(group_data.iloc[row_idx]['start']))
                    if key in variant_data:
                        sample_genotypes = variant_data[key]
                        merged_genotypes = combine_genotypes(merged_genotypes, sample_genotypes)
                    else:
                        print(f"Warning: Key {key} not found in variant_data.")

            # If the current kmer has exactly one match in the comparison matrix, skip adding this variant
            if [comparison[idx] for comparison in comparison_matrix].count(1) == 1:

                continue  # Skip this variant, don't add it to the MAP or PED data

            # Format genotype data
            for sample in samples:
                def format_genotype(gt):
                    if None in gt or -1 in gt:
                        return '0 0'
                    return ' '.join(map(str, [gt[0] + 1, gt[1] + 1]))
                genotype_str = format_genotype(merged_genotypes[sample])
                ped_data[sample].append(genotype_str)

            # Append to MAP data
            map_row = [chromosome, f"chr{chromosome}_grp{group}_idx{idx}", 0, group_data.iloc[0]['start']]
            map_data.append(map_row)

    # Write PED file
    ped_output_path = f"{output_prefix}.ped"
    with open(ped_output_path, 'w') as ped_file:
        for row in ped_data.values():
            ped_file.write('\t'.join(map(str, row)) + '\n')

    # Write MAP file
    map_output_path = f"{output_prefix}.map"
    with open(map_output_path, 'w') as map_file:
        for row in map_data:
            map_file.write('\t'.join(map(str, row)) + '\n')


def convert_to_plink_with_variants(config: Config):
    """Main function: replaces variant names and generates PLINK files."""
    csv_data = load_csv(config.grouped_variants_file)
    print("Initial Columns:", csv_data.columns.tolist())

    
    csv_data = process_comparison_column(csv_data)
    print("Columns after processing 'comparison':", csv_data.columns.tolist())

    
    if 'chromosome_group' not in csv_data.columns:
        raise ValueError("Missing 'chromosome_group' column in the input CSV.")

    
    variants = parse_fasta(config.ref_fasta)
    csv_data = replace_seq_with_variants(csv_data, variants)

    
    csv_data = extract_group_from_chromosome_group(csv_data)

    
    csv_data = add_start_column_if_missing(csv_data)

    
    if 'chromosome' not in csv_data.columns or 'position' not in csv_data.columns:
        def extract_chrom_pos(seq_id):
            match = re.match(r'Variant_(\d+)_(\d+)_(\d+)_(\d+)', seq_id)
            if match:
                chrom, _, pos, _ = match.groups()
                return chrom, pos
            else:
                return None, None
        csv_data[['chromosome', 'position']] = csv_data['sequence_id'].apply(
            lambda x: pd.Series(extract_chrom_pos(x))
        )

   
    updated_csv_path = os.path.join(config.output_dir, "updated_processed_comparison_results.csv")
    csv_data.to_csv(updated_csv_path, index=False)
    print(f"Updated CSV saved to {updated_csv_path}")

    
    output_prefix = os.path.join(config.output_dir, "plink_output")
    create_ped_and_map_files(csv_data, config.vcf_file, output_prefix)

    
    subprocess.run([
        'plink',
        '--file', output_prefix,
        '--make-bed',
        '--out', output_prefix
    ], check=True)

    print(f"PLINK files saved in {config.output_dir}")
