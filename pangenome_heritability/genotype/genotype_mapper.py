import os
import click
import pandas as pd
import re
import subprocess
from typing import Dict, NamedTuple
from config import Config  # Ensure this module contains the Config class

class PlinkFiles(NamedTuple):
    bed: str
    bim: str
    fam: str

def load_csv(file_path: str) -> pd.DataFrame:
    """Load the CSV file."""
    return pd.read_csv(file_path)

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
            # We ignore sequence lines as we only need variant names
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
        else:
            # If 'sequence_id' is already a variant name, leave it as is
            pass
    return csv_data

def process_comparison_column(csv_data: pd.DataFrame) -> pd.DataFrame:
    """Process the 'comparison' column to extract genotype information."""
    def process_comparison(s):
        try:
            comparison_list = eval(s)
            if len(comparison_list) == 1 and comparison_list[0] == 0:
                return [1]
            return comparison_list
        except Exception as e:
            print(f"Error processing comparison data: {s}. Error: {e}")
            return None
    csv_data['comparison'] = csv_data['comparison'].apply(process_comparison)
    csv_data = csv_data[csv_data['comparison'].notnull()]
    return csv_data

def create_ped_and_map_files(csv_data: pd.DataFrame, output_prefix: str):
    """Create PED and MAP files from the processed CSV data."""
    # Get list of samples
    samples = csv_data['sample'].unique()
    
    # Initialize PED data
    ped_data = {sample: [
        sample,  # Family ID
        sample,  # Individual ID
        '0',     # Paternal ID
        '0',     # Maternal ID
        '0',     # Sex (0=unknown)
        '-9'     # Phenotype (-9=missing)
    ] for sample in samples}
    
    # Initialize MAP data
    map_data = []

    # We need to collect genotype data per variant
    variants = csv_data['sequence_id'].unique()
    
    for idx, variant in enumerate(variants):
        map_row = [
            csv_data[csv_data['sequence_id'] == variant]['chromosome'].iloc[0],  # Chromosome
            variant,        # SNP ID
            '0',            # Genetic distance (can be zero)
            csv_data[csv_data['sequence_id'] == variant]['position'].iloc[0]  # Base-pair position
        ]
        map_data.append(map_row)
        
        # Now collect genotype data for each sample for this variant
        for sample in samples:
            sample_rows = csv_data[(csv_data['sample'] == sample) & (csv_data['sequence_id'] == variant)]
            if not sample_rows.empty:
                genotype_values = sample_rows['comparison'].iloc[0]
                # Assuming genotype_values is a list of integers
                # For simplicity, we'll encode '1' as '1 1' (homozygous), '0' as '0 0' (missing)
                genotype_str = ' '.join(['1' if val == 1 else '0' for val in genotype_values])
            else:
                genotype_str = '0 0'  # Missing genotype
            ped_data[sample].append(genotype_str)
    
    # Write PED file
    ped_output_path = f"{output_prefix}.ped"
    with open(ped_output_path, 'w') as f:
        for sample in samples:
            row = ped_data[sample]
            f.write('\t'.join(map(str, row)) + '\n')

    # Write MAP file
    map_output_path = f"{output_prefix}.map"
    with open(map_output_path, 'w') as f:
        for row in map_data:
            f.write('\t'.join(map(str, row)) + '\n')

def convert_to_plink_with_variants(config: Config):
    """Main function: replaces variant names and generates PLINK files."""
    # Load CSV
    csv_data = load_csv(config.grouped_variants_file)

    # Process 'comparison' column
    csv_data = process_comparison_column(csv_data)
    
    # Parse FASTA and replace 'sequence_id' with variant names
    variants = parse_fasta(config.ref_fasta)
    csv_data = replace_seq_with_variants(csv_data, variants)
    
    # Ensure 'chromosome' and 'position' columns are present
    if 'chromosome' not in csv_data.columns or 'position' not in csv_data.columns:
        # Extract from 'sequence_id' if possible
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
    
    # Save the updated CSV
    updated_csv_path = os.path.join(config.output_dir, "updated_processed_comparison_results.csv")
    csv_data.to_csv(updated_csv_path, index=False)
    print(f"Updated CSV saved to {updated_csv_path}")

    # Create PLINK PED and MAP files
    output_prefix = os.path.join(config.output_dir, "plink_output")

    create_ped_and_map_files(csv_data, output_prefix)

    # Run PLINK conversion
    subprocess.run([
        'plink',
        '--file', output_prefix,
        '--make-bed',
        '--out', output_prefix
    ], check=True)

    print(f"PLINK files saved in {config.output_dir}")