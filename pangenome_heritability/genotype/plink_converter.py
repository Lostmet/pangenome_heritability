import os
import subprocess
from typing import Dict, NamedTuple
import pandas as pd
from ..config import Config

class PlinkFiles(NamedTuple):
    bed: str
    bim: str
    fam: str

def convert_to_plink(config: Config, 
                    kmer_results: pd.DataFrame) -> PlinkFiles:
    """Convert k-mer comparison results to PLINK binary format"""
    output_prefix = os.path.join(config.output_dir, "variants")
    
    # Create PED file
    create_ped_file(kmer_results, f"{output_prefix}.ped")
    
    # Create MAP file
    create_map_file(kmer_results, f"{output_prefix}.map")
    
    # Convert to binary format
    subprocess.run([
        'plink',
        '--file', output_prefix,
        '--make-bed',
        '--out', output_prefix
    ], check=True)
    
    return PlinkFiles(
        bed=f"{output_prefix}.bed",
        bim=f"{output_prefix}.bim",
        fam=f"{output_prefix}.fam"
    )

def create_ped_file(kmer_results: pd.DataFrame, output_path: str):
    """Create PED file from k-mer results"""
    samples = kmer_results['sample'].unique()
    with open(output_path, 'w') as f:
        for sample in samples:
            sample_data = kmer_results[kmer_results['sample'] == sample]
            row = [
                sample,  # Family ID
                sample,  # Individual ID
                '0',    # Paternal ID
                '0',    # Maternal ID
                '0',    # Sex
                '-9'    # Phenotype
            ]
            row.extend(sample_data['genotype'].values)
            f.write('\t'.join(map(str, row)) + '\n')

def create_map_file(kmer_results: pd.DataFrame, output_path: str):
    """Create MAP file from k-mer results"""
    variants = kmer_results[['chrom', 'pos', 'window_id']].drop_duplicates()
    with open(output_path, 'w') as f:
        for _, var in variants.iterrows():
            row = [
                var['chrom'],
                var['window_id'],
                '0',
                str(var['pos'])
            ]
            f.write('\t'.join(row) + '\n')