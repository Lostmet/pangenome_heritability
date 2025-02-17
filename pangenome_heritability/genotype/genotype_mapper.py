import os
import click
import ast
import pandas as pd
import re
import subprocess
import pysam
from typing import Dict, NamedTuple
from ..config import Config  # Ensure this module contains the Config class
import numpy as np
import json
import glob


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
        #'SV': f"{output_prefix}.sv.vcf"#,
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
    #for variant_type, output_vcf in output_files.items():
    output_vcf = output_files['rSV']
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
                variant_type = 'rSV'
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


##############


# -----------------------------------
# Step 1: 处理 diff_array，生成 D_matrix
# -----------------------------------

def process_diff_array(input_csv: str, output_dir: str):
    """ 解析 diff_array 并提取 meta_array 的 pos，生成 D_matrix """

    input_path = os.path.abspath(input_csv)  # 将传入的相对路径转换为绝对路径

    if not os.path.exists(input_path):
        print(f"File not found: {input_path}")
        return

    df = pd.read_csv(input_path)

    group_dict = {}
    #group_meta_pos = {}  # 记录每个 group 的 seq1 的 meta_array 第一个 pos 值

    for _, row in df.iterrows():
        group_name = row['chromosome_group']
        seq_id = row['sequence_id']
        diff_array_str = row['diff_array']
        #meta_array_str = row['meta_array']  # 读取 meta_array

        try:
            diff_array = json.loads(diff_array_str)
        except json.JSONDecodeError:
            print(f"Error decoding diff_array for {group_name}. Skipping this row.")
            continue
        '''       
        try:
            meta_array = json.loads(meta_array_str.replace("'", "\""))
            if isinstance(meta_array, list) and len(meta_array) > 0:
                meta_pos_value = meta_array[0].get('pos', 'unknown')
            else:
                meta_pos_value = 'unknown'
        except json.JSONDecodeError:
            meta_pos_value = 'unknown'
        ''' 
        #if group_name not in group_meta_pos:
            #group_meta_pos[group_name] = meta_pos_value

        if diff_array != [1]:
            if group_name not in group_dict:
                group_dict[group_name] = []
            group_dict[group_name].append((seq_id, diff_array))

    for group_name, data in group_dict.items():
        seq_ids = [item[0] for item in data]
        diff_arrays = [item[1] for item in data]

        seq_count = len(diff_arrays[0])
        columns = ['sv_id'] + [f'rSV{i+1}' for i in range(seq_count)]

        rows = [[seq_id] + diff_array for seq_id, diff_array in zip(seq_ids, diff_arrays)]
        df_new = pd.DataFrame(rows, columns=columns)

        #meta_pos_value = group_meta_pos.get(group_name, 'unknown')
        output_file = os.path.join(output_dir, f"{group_name}_D_matrix.csv")
        df_new.to_csv(output_file, index=False)
        #print(f"Matrix for {group_name} saved to {output_file}")

# -----------------------------------
# Step 2: 读取 VCF，生成 X_matrix
# -----------------------------------

#2/8改了一下X矩阵，应该没啥问题）
def process_vcf_to_x_matrix(vcf_name: str, output_dir: str):
    """ 从 VCF 文件提取 GT 矩阵，并匹配 D_matrix 生成 X_matrix """
    vcf_file = pysam.VariantFile(vcf_name, 'r')  # 使用 pysam 打开 VCF 文件
    sample_names = list(vcf_file.header.samples)  # 获取样本名称列表

    # 存储 GT 数据的列表
    gt_data = []

    # 逐行遍历 VCF 文件
    for record in vcf_file:
        # 获取每个样本的 GT 信息
        gt_row = [record.contig, record.pos]  # 每行包含染色体和位置

        for sample in sample_names:
            gt = record.samples[sample]["GT"]
            #print(f"gt:{gt}")
            if gt == (None, None):  # 处理缺失值
                gt_row.append("./.")
            else:
                gt_row.append(f"{gt[0]}/{gt[1]}")  # 格式化为 '0/1'

        gt_data.append(gt_row)

    # 将 GT 数据转换为 DataFrame
    header = ["#CHROM", "POS"] + sample_names  # 只保留染色体、位置和样本名称
    df_vcf = pd.DataFrame(gt_data, columns=header)

    # 转换基因型函数
    def transform_gt(gt):
        if gt == './.': return 0
        elif gt == '0/0': return 0
        elif gt == '1/0' or gt == '0/1': return 1
        elif gt == '1/1': return 2
        else: return -1

    # 对样本列应用基因型转换
    for sample in sample_names:
        df_vcf[sample] = df_vcf[sample].apply(transform_gt)

    # 确保输出目录是绝对路径
    output_dir = os.path.abspath(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)  # 如果目录不存在，创建它

    # 从 output_dir 中读取 D_matrix 文件
    csv_files = [f for f in os.listdir(output_dir) if re.match(r'Group_\d+_\d+_\d+_D_matrix\.csv', f)]

    for csv_file in csv_files:
        match = re.match(r'Group_(\d+)_(\d+)_(\d+)_D_matrix\.csv', csv_file)
        if not match:
            continue

        chrom, number, pos = match.groups()
        pos = int(pos)
        chrom = str(chrom)

        # 过滤出匹配 CHROM 和 POS 的行
        vcf_row_index = df_vcf[(df_vcf["#CHROM"].astype(str) == chrom) & (df_vcf["POS"] == pos)].index

        if vcf_row_index.empty:
            print(f"Warning: #CHROM {chrom}, POS {pos} not found in VCF for {csv_file}")
            continue

        start_idx = vcf_row_index[0]
        csv_data = pd.read_csv(os.path.join(output_dir, csv_file), usecols=[0])
        num_rows_to_extract = len(csv_data)

        # 获取 VCF 子集数据
        vcf_subset = df_vcf.iloc[start_idx : start_idx + num_rows_to_extract].reset_index(drop=True)
        gt_matrix = vcf_subset[sample_names]

        # 合并 CSV 数据和 GT 矩阵
        csv_data.columns = [f"{chrom}_{pos}"]
        csv_data = pd.concat([csv_data, gt_matrix], axis=1)

        # 确保 X_matrix 输出到指定目录
        output_file = os.path.join(output_dir, csv_file.replace("_D_matrix.csv", "_X_matrix.csv"))
        csv_data.to_csv(output_file, index=False)

        if not os.path.exists(output_file):
            print(f"Warning: {output_file} not created!")
        else:
            continue
            #print(f"Saved modified CSV: {output_file}")
    return sample_names
# -----------------------------------
# Step 3: 计算 T_matrix
# -----------------------------------

def compute_t_matrix(output_dir: str):
    """ 计算 D_matrix × X_matrix 并保存 T_matrix """
    
    # 确保输出目录是绝对路径
    output_dir = os.path.abspath(output_dir)

    # 获取 output_dir 下的所有文件
    files = os.listdir(output_dir)
    d_files = [f for f in files if re.match(r'Group_\d+_\d+_\d+_D_matrix\.csv', f)]
    x_files = [f for f in files if re.match(r'Group_\d+_\d+_\d+_X_matrix\.csv', f)]

    # 排序文件
    d_files.sort()
    x_files.sort()

    # 确保 D_matrix 和 X_matrix 文件一一对应
    if len(d_files) != len(x_files):
        print("Warning: D_matrix and X_matrix file counts do not match.")
        return

    for d_file, x_file in zip(d_files, x_files):
        match_d = re.match(r'Group_(\d+)_(\d+)_(\d+)_D_matrix\.csv', d_file)
        match_x = re.match(r'Group_(\d+)_(\d+)_(\d+)_X_matrix\.csv', x_file)

        if not match_d or not match_x:
            print(f"Skipping unmatched files: {d_file}, {x_file}")
            continue

        # 使用完整的路径读取文件
        df_d = pd.read_csv(os.path.join(output_dir, d_file))
        df_x = pd.read_csv(os.path.join(output_dir, x_file))

        d_data = df_d.iloc[:, 1:].values
        x_data = df_x.iloc[:, 1:].values

        # 计算 T_matrix
        t_matrix = np.dot(d_data.T, x_data)

        # 创建 T 矩阵并转置
        t_df = pd.DataFrame(t_matrix, index=df_d.columns[1:], columns=df_x.columns[1:])

        # 获取 X_matrix 的第一列名字，并赋值给 T_matrix 的 index name
        first_column_name = df_x.columns[0]
        t_df.index.name = first_column_name

        # 保存 T_matrix 到指定目录
        t_matrix_file = os.path.join(output_dir, d_file.replace("_D_matrix.csv", "_T_matrix.csv"))
        t_df.to_csv(t_matrix_file, index=True, header=True)
        #print(f"Saved T matrix: {t_matrix_file}")

def save_rsv_meta(input: str, output:str):
    # 读取原始CSV文件
    df = pd.read_csv(input)

    # 找到第一列中重复的项，并保留所有重复项
    df = df[df.duplicated(subset=df.columns[0], keep=False)]

    #print("删除唯一项并保留重复项成功！")

    # 定义一个集合来存储所有提取的meta_array数据（使用集合去除重复项）
    extracted_meta_data = set()
    # **只遍历包含 'seq1' 的行**
    df_filtered = df[df[df.columns[1]] == 'seq1']  # 过滤出第二列为 'seq1' 的行
    # 遍历所有行
    for index, row in df_filtered.iterrows():
        group_name = row[df_filtered.columns[0]]
        
        # 解析 diff_array 和 meta_array
        diff_array = ast.literal_eval(row['diff_array'])
        meta_array = ast.literal_eval(row['meta_array'])
        
        # 遍历 diff_array，找到值为 1 的位置，提取对应的 meta_array
        for i in range(len(diff_array)):
            extracted_meta_data.add((group_name, str(meta_array[i])))

    # 将提取的 meta_array 数据转换为列表
    meta_list = [(group, ast.literal_eval(item)) for group, item in extracted_meta_data]

    # 按照 pos 字段排序
    meta_list.sort(key=lambda x: x[1]['pos'])

    # 创建 DataFrame
    final_data = pd.DataFrame(meta_list, columns=["group_name", "meta_array"])

    # 去重处理
    final_data['meta_array_str'] = final_data['meta_array'].apply(lambda x: str(x))
    final_data.drop_duplicates(subset=['meta_array_str'], inplace=True)

    # 删除辅助列
    final_data.drop(columns=['meta_array_str'], inplace=True)

    # **直接在内存中进行后缀添加**
    final_data['group_name_suffix'] = final_data.groupby('group_name').cumcount() + 1
    final_data['group_name'] = final_data['group_name'] + '_rSV' + final_data['group_name_suffix'].astype(str)

    # 删除临时列
    final_data.drop(columns=['group_name_suffix'], inplace=True)

    # 保存文件到指定目录（固定文件名为 rsv_meta.csv）
    file_name = "rsv_meta.csv"
    output_file = os.path.join(output, file_name)
    final_data.to_csv(output_file, index=False)

    print(f"{file_name} generated！")

######下面提取T矩阵对应的GT矩阵

def find_t_matrix_file(chrom, number, out):
    """ 查找符合 Group_{chrom}_{number}_*_T_matrix.csv 形式的文件 """
    pattern = os.path.join(out, f"Group_{chrom}_{number}_*_T_matrix.csv")  # 拼接路径和模式
    matching_files = glob.glob(pattern)  # 匹配文件
    return matching_files[0] or None  # 如果没有找到文件，返回空列表
        # 确保输出目录是绝对路径
    #output_dir = os.path.abspath(output_dir)

    # 获取 output_dir 下的所有文件
    #files = os.listdir(output_dir)
    #d_files = [f for f in files if re.match(r'Group_\d+_\d+_\d+_D_matrix\.csv', f)]
    # 如果当前目录找不到，递归搜索子目录
    #if not matching_files:
        #matching_files = glob.glob(f"**/{pattern}", recursive=True)

    

def convert_gt(value):
    """ 将数值转换为 VCF 格式的 GT 类型 """
    try:
        value = int(value)
        if value == 0:
            return "0/0"
        elif value == 1:
            return "1/0"
        elif value == 2:
            return "1/1"
        else:
            return "./."
    except:
        return "./."  # 遇到异常情况时填充缺失数据

def extract_vcf_sample(input_csv: str, output_vcf: str, out: str):
    """ 从 rsv_meta.csv 解析 group_name，并提取 GT 样本数据 """
    #input_csv = "rsv_meta_暴力提取.csv"
    #output_vcf = "vcf_samples.vcf"

    df_meta = pd.read_csv(input_csv)

    gt_data = []

    for index, row in df_meta.iterrows():
        group_name = row["group_name"]  # 例："Group_2_19_某pos_rSV1"
        match = re.match(r"Group_(\d+)_(\d+)_(\d+)_([rR][sS][vV]\d+)", group_name)

        if not match:
            print(f"⚠️ Skipping invalid group_name: {group_name}")
            continue

        chrom, number, pos, rsv_name = match.groups()

        # **查找匹配的 T_matrix 文件**
        t_matrix_file = find_t_matrix_file(chrom, number, out)

        if not t_matrix_file:
            print(f"⚠️ Warning: No matching T_matrix found for {group_name} (expected: Group_{chrom}_{number}_*_T_matrix.csv)")
            continue

        #print(f"✅ Found T_matrix file: {t_matrix_file} for {group_name}")

        # **读取 T_matrix**
        df_t = pd.read_csv(t_matrix_file, index_col=0)

        if rsv_name not in df_t.index:
            print(f"⚠️ Warning: {rsv_name} not found in {t_matrix_file}")
            continue

        # **提取 GT 行，并转换为 VCF GT 格式**
        gt_values = df_t.loc[rsv_name].apply(convert_gt).tolist()
        gt_data.append(gt_values)

    # **转换为 DataFrame**
    gt_df = pd.DataFrame(gt_data)

    # **保存 VCF（无表头）**
    gt_df.to_csv(output_vcf, index=False, header=False, sep="\t")

    #print(f"✅ VCF GT 数据已保存至 {output_vcf}！")

########生成rsv

def vcf_generate(sample_names: list, csv_file: str, output_vcf: str, gt_file: str, output_filled_vcf: str):
    # 读取 CSV 文件
    df = pd.read_csv(csv_file)

    # 解析 group_name 以提取 #CHROM 和编号
    vcf_data = []

    for index, row in df.iterrows():
        try:
            # 解析 group_name
            group_name = row["group_name"]  # ID 直接设为 group_name
            group_name_parts = group_name.split("_")
            chrom = group_name_parts[1]  # 提取 CHROM
            
            # 解析 meta_array 字段
            meta_data = ast.literal_eval(row["meta_array"])
            pos = meta_data["pos"]  # 解析 POS
            ref = meta_data["ref"] if meta_data["ref"] else "_" # 解析 REF，空为"_"
            alt = meta_data["alt"] if meta_data["alt"] else "_"  # 解析 ALT，若为空则设为 "_"

            # 组装 VCF 结构
            vcf_data.append([chrom, pos, group_name, ref, alt, ".", ".", "TYPE=rSV", "GT"])
        except Exception as e:
            print(f"Error processing row {index}: {e}")
            continue

    # 转换为 DataFrame
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    vcf_df = pd.DataFrame(vcf_data, columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"])

    # **按 CHROM 和 POS 进行排序**
    vcf_df.sort_values(by=["#CHROM", "POS"], inplace=True)
    #print(vcf_df.columns)
    """ 从原始 VCF 提取列名，并保持一致写入新 VCF """
    with open(output_vcf, "w") as outfile:
        header = "\t".join(columns) + "\t" + "\t".join(sample_names)  # 提取正确的列名
        outfile.write("##fileformat=VCFv4.2\n")  # VCF 版本信息
        outfile.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        outfile.write(header + "\n")  # 写入正确的表头

        vcf_df.to_csv(outfile, sep="\t", index=False, header=False)
    #print(f"✅ VCF 头部提取完成，已写入 {output_vcf}")

    # 读取 VCF 和 GT 数据
    with open(output_vcf, "r") as vcf, open(gt_file, "r") as gt, open(output_filled_vcf, "w") as out_vcf:
        for line in vcf:
            if line.startswith("#"):  # 头部信息保持不变
                out_vcf.write(line)
            else:
                gt_line = gt.readline().strip()  # 读取 GT 数据行
                fields = line.strip().split("\t")  # 解析 VCF 行
                fields[9:] = gt_line.split()  # **从第 10 列 (索引 9) 开始替换 GT**
                out_vcf.write("\t".join(fields) + "\n")
    #print(f"✅ GT 数据填充完成，最终 VCF 文件已保存至 {output_filled_vcf}！")

####检测错误GT
def detect_abnormal(out: str):
    # 获取当前目录下所有 T_matrix 文件
    t_matrix_files = [f for f in os.listdir(out) if re.match(r'Group_\d+_\d+_\d+_T_matrix\.csv', f)]
    # 统计信息列表
    stats_list = []
    # 添加前缀 'out/' 到每个文件名
    directory = out  # 'out' 是你要加前缀的目录
    t_matrix_files = [os.path.join(directory, file_name) for file_name in t_matrix_files]
    #print(t_matrix_files)
    # 遍历所有 T_matrix 文件
    for t_file in t_matrix_files:
        df = pd.read_csv(t_file, index_col=0)  # 读取 T_matrix，第一列作为索引

        # 计算 >0 的数量
        num_gt_0 = (df.values > 0).sum()

        # 计算 >2 的数量
        num_gt_2 = (df.values > 2).sum()

        # 计算异常率
        abnormal_rate = num_gt_2 / num_gt_0 if num_gt_0 > 0 else 0

        # 记录统计结果
        stats_list.append([t_file, num_gt_0, num_gt_2, abnormal_rate])

    # 转换为 DataFrame 并保存 T_matrix_abnormal.csv
    df_stats = pd.DataFrame(stats_list, columns=["group_name", "number > 0", "number >2", "abnormal_rate"])
    output_file = os.path.join(out, "T_matrix_abnormal.csv")
    df_stats.to_csv(output_file, index=False)
    print("Saved: T_matrix_abnormal.csv")

    # 计算所有 T_matrix 的汇总统计
    total_gt_0 = df_stats["number > 0"].sum()
    total_gt_2 = df_stats["number >2"].sum()
    total_abnormal_rate = total_gt_2 / total_gt_0 if total_gt_0 > 0 else 0

    # 生成 T_matrix_abnormal_all.csv
    df_total = pd.DataFrame([[total_gt_0, total_gt_2, total_abnormal_rate]], 
                            columns=["number > 0", "number >2", "abnormal_rate"])
    output_file = os.path.join(out, "T_matrix_abnormal_all.csv")

    df_total.to_csv(output_file, index=False)
    print("Saved: T_matrix_abnormal_all.csv")


