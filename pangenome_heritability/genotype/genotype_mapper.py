import os
import click
import pandas as pd
import pysam
import subprocess
from typing import Dict, NamedTuple
from config import Config  # 从外部模块导入 Config 类


class PlinkFiles(NamedTuple):
    bed: str
    bim: str
    fam: str


def load_csv(file_path: str) -> pd.DataFrame:
    """加载 CSV 文件"""
    return pd.read_csv(file_path)


def parse_fasta(file_path: str) -> Dict[str, list]:
    """解析 FASTA 文件，提取变异名及其对应的序列"""
    variants = {}
    current_group = None
    current_variant = None
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:]
                if "Variant" in header:
                    current_variant = header
                else:
                    current_group = header
                    if current_group not in variants:
                        variants[current_group] = []
            else:
                if current_group and current_variant:
                    variants[current_group].append(current_variant)
                    current_variant = None
    return variants


def replace_seq_with_variants(csv_data: pd.DataFrame, variants: Dict[str, list]) -> pd.DataFrame:
    """将 sequence_id 替换为变异名"""
    for index, row in csv_data.iterrows():
        group = row['chromosome_group']
        if 'seq' in row['sequence_id']:
            seq_id = int(row['sequence_id'].replace('seq', '')) - 1
            if group in variants and seq_id < len(variants[group]):
                csv_data.at[index, 'sequence_id'] = variants[group][seq_id]
            else:
                print(f"No variant available for {group} at sequence {seq_id}")
        else:
            print(f"Invalid sequence_id format: {row['sequence_id']}")
    return csv_data


def create_ped_file(kmer_results: pd.DataFrame, output_path: str):
    """创建 PED 文件"""
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
            row.extend(sample_data['genotype'].fillna('0 0').values)
            f.write('\t'.join(map(str, row)) + '\n')


def create_map_file(kmer_results: pd.DataFrame, output_path: str):
    """创建 MAP 文件"""
    variants = kmer_results[['chrom', 'pos', 'window_id']].drop_duplicates()
    with open(output_path, 'w') as f:
        for _, var in variants.iterrows():
            row = [
                var['chrom'],  # Chromosome
                var['window_id'],  # SNP ID
                '0',  # Genetic distance
                var['pos']  # Base-pair position
            ]
            f.write('\t'.join(map(str, row)) + '\n')


def convert_to_plink_with_variants(config: Config):
    """主流程: 替换变异名并生成 PLINK 文件"""
    # 加载 CSV 和 FASTA
    csv_data = load_csv(config.grouped_variants_file)
    variants = parse_fasta(config.ref_fasta)

    # 替换 sequence_id 为 variants
    updated_csv = replace_seq_with_variants(csv_data, variants)

    # 保存更新后的 CSV
    updated_csv_path = os.path.join(config.output_dir, "updated_processed_comparison_results.csv")
    updated_csv.to_csv(updated_csv_path, index=False)
    print(f"Updated CSV saved to {updated_csv_path}")

    # 创建 PLINK PED 和 MAP 文件
    output_prefix = os.path.join(config.output_dir, "plink_output")
    create_ped_file(updated_csv, f"{output_prefix}.ped")
    create_map_file(updated_csv, f"{output_prefix}.map")

    # 运行 PLINK 转换
    subprocess.run([
        'plink',
        '--file', output_prefix,
        '--make-bed',
        '--out', output_prefix
    ], check=True)

    print(f"PLINK files saved in {config.output_dir}")