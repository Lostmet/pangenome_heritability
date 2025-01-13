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
    """Process diff_array column to extract genotype information"""
    print("Columns before processing 'diff_array':", csv_data.columns.tolist())
    
    if 'diff_array' not in csv_data.columns:
        print("Warning: 'diff_array' column is missing or empty.")
        return csv_data

    def process_diff_array(s):
        try:
            cleaned_s = re.sub(r',\s*$', '', s)
            return eval(cleaned_s)
        except Exception as e:
            print(f"Error processing diff_array data: {s}. Error: {e}")
            return None

    csv_data['diff_array'] = csv_data['diff_array'].apply(process_diff_array)
    csv_data = csv_data[csv_data['diff_array'].notnull()]
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
    
    
    variant_data = {}
    for record in vcf:
        chromosome = str(record.contig)
        position = str(record.pos)
        key = (chromosome, position)
        variant_data[key] = {sample: record.samples[sample]['GT'] for sample in samples}
    vcf.close()

    
    def combine_genotypes(existing, new):
        def combine_gt(gt1, gt2):
            if None in gt1 or -1 in gt1:
                return gt2
            if None in gt2 or -1 in gt2:
                return gt1
            return tuple(max(a, b) for a, b in zip(gt1, gt2))
        return {sample: combine_gt(existing.get(sample, (0, 0)), new.get(sample, (0, 0))) 
                for sample in samples}

    
    ped_data = {sample: [sample, sample, '0', '0', '0', '-9'] for sample in samples}
    map_data = []

    
    for _, row in csv_data.iterrows():
        try:
            
            chromosome = row['chromosome_group'].split('_')[1]
            group = row['chromosome_group'].split('_')[-1]
            
            
            meta_array = eval(row['meta_array'])
            diff_array = eval(row['diff_array'])
            
            
            for i, (meta, diff) in enumerate(zip(meta_array, diff_array)):
                pos = str(meta['pos'])
                ref = meta['ref']
                alt = meta['alt']
                
                
                vcf_key = (str(chromosome), pos)
                if vcf_key not in variant_data:
                    continue

                
                variant_count = diff_array.count(1)
                variant_type = "SV" if variant_count == 1 else "RSV"
                
                
                variant_id = f"{variant_type}_chr{chromosome}_grp{group}_pos{pos}_{ref}_{alt}"
                
               
                map_row = [chromosome, variant_id, '0', pos]
                map_data.append(map_row)
                
                
                for sample in samples:
                    gt = variant_data[vcf_key].get(sample, (0, 0))
                   
                    if None in gt or -1 in gt:
                        genotype_str = '0 0'
                    else:
                        genotype_str = ' '.join(map(str, [gt[0] + 1, gt[1] + 1]))
                    ped_data[sample].append(genotype_str)
                    
        except Exception as e:
            print(f"Error processing row: {str(e)}")
            continue

    
    with open(f"{output_prefix}.ped", 'w') as f:
        for sample, data in ped_data.items():
            f.write('\t'.join(map(str, data)) + '\n')

    
    with open(f"{output_prefix}.map", 'w') as f:
        for row in map_data:
            f.write('\t'.join(map(str, row)) + '\n')

def create_vcf_file(csv_data: pd.DataFrame, vcf_file: str, grouped_variants: str, output_prefix: str):
    """Generate new VCF file based on CSV data and original VCF"""
    import pysam
    
    # Open VCF file and get sample information
    vcf = pysam.VariantFile(vcf_file)
    samples = list(vcf.header.samples)
    
    # Create output files for SV and rSV
    output_files = {
        'SV': f"{output_prefix}.sv.vcf",
        'rSV': f"{output_prefix}.rsv.vcf"
    }
    
    # Extract variant data from VCF
    variant_data = {}
    print("Reading VCF file...")
    for record in vcf:
        chromosome = str(record.contig)
        position = int(record.pos)
        key = (chromosome, position)
        variant_data[key] = record
    vcf.close()
    
    # Create files for each variant type and write headers
    for variant_type, output_vcf in output_files.items():
        with open(output_vcf, 'w') as f:
            f.write('##fileformat=VCFv4.2\n')
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples) + '\n')
    
    # Process data by group
    for group_name, group_data in csv_data.groupby('chromosome_group'):
        try:
            #print(f"Processing group: {group_name}")
            
            # Get meta_array and diff_array
            if 'meta_array' not in group_data.columns or 'diff_array' not in group_data.columns:
                #print(f"Group {group_name} missing required columns")
                continue
                
            meta_array = eval(group_data.iloc[0]['meta_array'])
            diff_array_matrix = [eval(row['diff_array']) for _, row in group_data.iterrows()]
            
            # Get variant information
            variants_info = get_variants_info(grouped_variants, group_name)
            if not variants_info:
                continue
            
            # Process variants for each column
            for col_idx in range(len(meta_array)):
                # Get diff values for this column across all rows
                col_diffs = [diff[col_idx] if col_idx < len(diff) else 0 for diff in diff_array_matrix]
                variant_count = sum(col_diffs)
                
                # Determine variant type based on variant_count
                variant_type = 'SV' if variant_count == 1 else 'rSV'
                output_vcf = output_files[variant_type]
                
                # Get variant information
                meta = meta_array[col_idx]
                ref = meta['ref']
                alt = meta['alt']
                
                # Get variants marked as 1
                variant_indices = [i for i, diff in enumerate(col_diffs) if diff == 1]
                if not variant_indices:
                    continue
                
                # Use first variant position as representative
                chrom, pos = variants_info[variant_indices[0]]
                key = (str(chrom), int(pos))
                if key not in variant_data:
                    continue
                
                record = variant_data[key]
                
                # Format genotype data, replace None with "."
                formatted_genotypes = []
                for sample in samples:
                    if sample in record.samples:
                        gt = record.samples[sample]['GT']
                        if gt is None:
                            formatted_genotypes.append('./.')
                        else:
                            formatted_genotypes.append('/'.join(map(str, [x if x is not None else '.' for x in gt])))
                    else:
                        formatted_genotypes.append('./.')
                
                # Write VCF line
                with open(output_vcf, 'a') as f:
                    vcf_line = [
                        chrom,
                        str(pos),
                        f"{variant_type}_chr{chrom}_grp{group_name}_pos{pos}",
                        ref,
                        alt,
                        '.',
                        'PASS',
                        f"TYPE={variant_type};GROUP={group_name}",
                        'GT'
                    ] + formatted_genotypes
                    
                    f.write('\t'.join(map(str, vcf_line)) + '\n')
                
        except Exception as e:
            print(f"Error processing group {group_name}: {str(e)}")
            continue

    print("VCF file generation complete")

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

def parse_variant_id(variant_id):
    """Parse variant information from sequence ID in variants.fasta
    Example: Variant_2_2_50381_50381 -> (2, 50381)
    """
    parts = variant_id.split('_')
    if len(parts) >= 5 and parts[0] == 'Variant':
        chromosome = parts[1]
        pos = int(parts[3])
        return chromosome, pos
    return None, None

def get_variants_info(fasta_file, group_name):
    """Get variant information for specified group from variants.fasta
    Example:
    >Group_2_1
    ACGT...
    >Variant_2_1_50381_50381
    ACGT...
    >Variant_2_1_50382_50382
    ACGT...
    """
    variants = []
    group_found = False
    
    #print(f"Reading variant information for {group_name} from {fasta_file}")
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                header = line.strip()[1:]
                #print(f"Reading header: {header}")
                
                if header == group_name:
                    group_found = True
                    #print(f"Found target group: {group_name}")
                    continue
                
                if group_found:
                    if header.startswith('Group_'):
                        #print(f"Encountered next group, ending current group processing")
                        break
                    elif header.startswith('Variant_'):
                        # Parse variant ID, example: Variant_2_1_50381_50381
                        parts = header.split('_')
                        if len(parts) >= 5:
                            chrom = parts[1]
                            pos = int(parts[3])
                            variants.append((chrom, pos))
                            #print(f"Added variant: chr{chrom}:{pos}")
    
    #print(f"Found {len(variants)} variants for group {group_name}")
    return variants
