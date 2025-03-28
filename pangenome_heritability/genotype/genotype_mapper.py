import os
import click
import ast
import pandas as pd
import re
import subprocess
import pysam
from typing import Dict, NamedTuple
from ..config import Config  
import numpy as np
import json
import glob
from concurrent.futures import ThreadPoolExecutor
from ..utils.logging_utils import get_logger, log_tqdm_summary


logger = get_logger(__name__)



class PlinkFiles(NamedTuple):
    bed: str
    bim: str
    fam: str

def load_csv(file_path: str) -> pd.DataFrame:
    """Load the CSV file with explicit handling for column names."""
    try:
        csv_data = pd.read_csv(file_path, sep=',', encoding='utf-8')
    except Exception as e:
        logger.warning(f"Error reading CSV file: {e}")
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
                logger.warning(f"No variant available for {group} at sequence {seq_id}")
    return csv_data

def process_comparison_column(csv_data: pd.DataFrame) -> pd.DataFrame:
    """Process diff_array column to extract genotype information"""
    click.echo("Columns before processing 'diff_array':", csv_data.columns.tolist())
    
    if 'diff_array' not in csv_data.columns:
        click.echo("Warning: 'diff_array' column is missing or empty.")
        return csv_data

    def process_diff_array(s):
        try:
            cleaned_s = re.sub(r',\s*$', '', s)
            return eval(cleaned_s)
        except Exception as e:
            click.echo(f"Error processing diff_array data: {s}. Error: {e}")
            return None

    csv_data['diff_array'] = csv_data['diff_array'].apply(process_diff_array)
    csv_data = csv_data[csv_data['diff_array'].notnull()]
    return csv_data

# -----------------------------------
# Step 1: 处理 diff_array，生成 D_matrix
# -----------------------------------
def process_diff_array(input_csv: str, output_dir: str):
    """ 解析 diff_array 并提取 meta_array 的 pos，生成 D_matrix """
    input_path = os.path.abspath(input_csv)  # 将传入的相对路径转换为绝对路径

    if not os.path.exists(input_path):
        logger.warning(f"File not found: {input_path}")
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
            logger.warning(f"Error decoding diff_array for {group_name}. Skipping this row.")
            continue

        if group_name not in group_dict:
            group_dict[group_name] = []
        group_dict[group_name].append((seq_id, diff_array[1:])) # 去掉diff_array第一项

    for group_name, data in group_dict.items():
        seq_ids = [item[0] for item in data]
        diff_arrays = [item[1] for item in data]

        seq_count = len(diff_arrays[0])
        columns = ['sv_id'] + [f'rSV{i+1}' for i in range(seq_count)]

        rows = [[seq_id] + diff_array for seq_id, diff_array in zip(seq_ids, diff_arrays)]
        df_new = pd.DataFrame(rows, columns=columns)

        output_file = os.path.join(output_dir, f"{group_name}_D_matrix.csv")
        df_new.to_csv(output_file, index=False)

# -----------------------------------
# Step 2: 读取 VCF，生成 X_matrix
# -----------------------------------
def sample_name_contract(vcf_name: str):
    """ 从 VCF 文件提取 GT 矩阵，并匹配 D_matrix 生成 X_matrix """
    vcf_file = pysam.VariantFile(vcf_name, 'r')  # 使用 pysam 打开 VCF 文件
    sample_names = list(vcf_file.header.samples)  # 获取样本名称列表
    return sample_names

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
        if gt == './.': return -999 # 改为缺失值
        elif gt == '0/0': return 0
        elif gt == '1/0' or gt == '0/1': return 1
        elif gt == '1/1': return 2
        else: return -1

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
            logger.warning(f"Warning: #CHROM {chrom}, POS {pos} not found in VCF for {csv_file}")
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
            logger.warning(f"Warning: {output_file} not created!")
        else:
            continue

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
        logger.warning("Warning: D_matrix and X_matrix file counts do not match.")
        return

    for d_file, x_file in zip(d_files, x_files):
        match_d = re.match(r'Group_(\d+)_(\d+)_(\d+)_D_matrix\.csv', d_file)
        match_x = re.match(r'Group_(\d+)_(\d+)_(\d+)_X_matrix\.csv', x_file)

        if not match_d or not match_x:
            logger.warning(f"Warning: Skipping unmatched files: {d_file}, {x_file}")
            continue

        # 使用完整的路径读取文件
        df_d = pd.read_csv(os.path.join(output_dir, d_file))
        df_x = pd.read_csv(os.path.join(output_dir, x_file))

        d_data = df_d.iloc[:, 1:].values
        x_data = df_x.iloc[:, 1:].values
        # 计算 T_matrix
        t_matrix = np.dot(d_data.T, x_data)

        # 创建 T 矩阵
        t_df = pd.DataFrame(t_matrix, index=df_d.columns[1:], columns=df_x.columns[1:])

        # 获取 X_matrix 的第一列名字，并赋值给 T_matrix 的 index name
        first_column_name = df_x.columns[0]
        t_df.index.name = first_column_name

        # 保存 T_matrix 到指定目录
        t_matrix_file = os.path.join(output_dir, d_file.replace("_D_matrix.csv", "_T_matrix.csv"))
        t_df.to_csv(t_matrix_file, index=True, header=True)

import concurrent.futures
from tqdm import tqdm

def process_row(row):
    # 将字符串转换为列表
    meta_list = ast.literal_eval(row['meta_array'])

    # 循环从第二个元素开始修改
    for j in range(1, len(meta_list)):
        # 获取前一个元素的ref的最后一个字符
        prev_ref_last_char = meta_list[j-1]['ref'][-1]
        
        # 更新当前元素的pos, ref, alt
        meta_list[j]['pos'] -= 1
        meta_list[j]['ref'] = prev_ref_last_char + meta_list[j]['ref']
        meta_list[j]['alt'] = prev_ref_last_char + meta_list[j]['alt']
    
    return meta_list

def save_rSV_meta(input: str, output: str, threads: int):
    # 读取原始CSV文件
    df = pd.read_csv(input)
    # 定义集合存储提取的元数据
    extracted_meta_data = set()
    click.echo("Pre-processing meta info of rSVs...")
    # 使用多线程并行化处理meta_array列，并加上进度条
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        # 使用tqdm包装executor.map，显示进度条
        results = list(tqdm(executor.map(process_row, [row for _, row in df.iterrows()]), total=len(df), desc="Processing info", unit="SV"))
    click.echo("Done")
    # 更新DataFrame中的meta_array列为修改后的结果
    for i, result in enumerate(results):
        df.at[i, 'meta_array'] = str(result)
    
    # 按group分组处理数据
    grouped = df.groupby(df.columns[0])

    for group_name, group_df in grouped:
        diff_arrays = []
        meta_arrays = []
        # 收集当前group所有seq的diff和meta
        for _, row in group_df.iterrows():
            try:
                diff = ast.literal_eval(row['diff_array'])
                meta = ast.literal_eval(row['meta_array'])
                diff_arrays.append(diff)
                meta_arrays.append(meta)
            except:
                continue
        
        # 确保存在有效数据并确定数组长度
        if len(diff_arrays) == 0:
            continue
        num_positions = len(diff_arrays[0])
        
        # 遍历每个位置（基因位点）
        for i in range(num_positions):
            # 检查该位置是否在任何seq中有diff=1
            for seq_idx in range(len(diff_arrays)):
                if diff_arrays[seq_idx][i] == 1:
                    # 选择第一个有效meta并去重
                    extracted_meta_data.add((group_name, str(meta_arrays[seq_idx][i])))
                    break  # 找到第一个有效的即可
    # 数据结构转换与预处理
    meta_list = [(group, ast.literal_eval(item)) for group, item in extracted_meta_data]
    click.echo("Sorting rSVs...")
    # 排序，先利用染色体，再用组别号，再用pos号
    meta_list.sort(key=lambda x: (int(x[0].split('_')[1]), int(x[0].split('_')[2]), x[1]['pos']))
    click.echo("Creating DataFrame...")
    # 创建DataFrame
    final_data = pd.DataFrame(meta_list, columns=["group_name", "meta_array"])
    click.echo("Grouping rSVs by group_name...")
    grouped_final_data = final_data.groupby(final_data.columns[0])
    click.echo("Done")
    # 定义处理函数，用于每个group_name的计算
    def rSV_info_modify(group_df):
        meta_arrays = []  # 存旧的
        ins = 0  # 这必须是我的神来之笔
        for _, row in group_df.iterrows():
            meta = row['meta_array']
            meta_arrays.append(meta)
        
        for i in range(len(meta_arrays)):  # 按照顺序处理refaltpos
            ref = meta_arrays[i]['ref']
            alt = meta_arrays[i]['alt']
            pos = meta_arrays[i]['pos']
            if ref[-1] == "-" and alt[0] != "-":
                pos -= ins
                ins += len(ref) - 1  # 给这个insertion对后面的进行清理
                ref = ref[0]  # ref的占位符去除
            elif alt[0] == "-" and alt[1] != "-": # poly_ins特有的'ref': '-----------', 'alt': '-TATATATATA'}"
                pos -= ins
                ins += len(ref) - 1
                ref = meta_arrays[i-1]['ref'][0] # ref=上一个的ref的第一位
                alt = meta_arrays[i-1]['ref'][0] + alt[1:]
            elif ref[0] == "-":  # 因为前面是insertion导致受损的变异
                pos -= ins
                ref = meta_arrays[i-1]['ref'][0] + ref[1:]
                alt = meta_arrays[i-1]['ref'][0]
            else:  # 不然就是deletion了，直接保留第一位就行
                pos -= ins
                alt = alt[0]
            meta_arrays[i]['ref'] = ref  # 倒反天罡替换一下
            meta_arrays[i]['alt'] = alt
            meta_arrays[i]['pos'] = pos

    # 使用ThreadPoolExecutor并行化group处理
    with ThreadPoolExecutor(max_workers=threads) as executor:
        # 提交任务并显示进度条
        list(tqdm(executor.map(rSV_info_modify, [group_df for _, group_df in grouped_final_data]),
                             desc="Converting to VCF v4.2 format...", unit="group",total=len(grouped_final_data)))

    #print(f'final_data: {final_data}')
    '''
    final_data:             group_name                                         meta_array
0      Group_1_1_17122  {'pos': 17122, 'ref': 'ATGTCACAGCCCGTCCCGGAATT...
1      Group_1_1_17122  {'pos': 17916, 'ref': 'T', 'alt': 'TTAACGGAAAC...
2      Group_1_1_17122  {'pos': 17916, 'ref': 'TTAACGGTACCGTTGGGTTTTTA...
3      Group_1_2_26725  {'pos': 26725, 'ref': 'CGCACCTTTAGTTATGTCACGCC...
4      Group_1_2_26725  {'pos': 31058, 'ref': 'A', 'alt': 'ATCAAGAACAA...
    '''
    # 添加后缀编号 # .groupby这里可以调整输入的是第一步还是第二步处理
    final_data['group_name_suffix'] = final_data.groupby('group_name').cumcount() + 1
    final_data['group_name'] = final_data['group_name'] + '_rSV' + final_data['group_name_suffix'].astype(str)
    final_data = final_data.drop(columns=['group_name_suffix'])

    # 保存文件
    file_name = "rSV_meta.csv"
    output_file = os.path.join(output, file_name)
    final_data.to_csv(output_file, index=False)
    click.echo(f"{file_name} generated!")


######下面提取T矩阵对应的GT矩阵
from concurrent.futures import ThreadPoolExecutor
def convert_gt(value):
    """ 将数值转换为 VCF 格式的 GT 类型 """
    try:
        value = int(value)
        if value == 0:
            return "0/0"
        elif value == 1:
            return "1/0"
        elif value >= 2:
            return "1/1" # 都当2
        else:
            return "./." 
    except:
        return "./."  # 遇到异常情况时填充缺失数据

import os
import re
import pandas as pd
import click
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict

def find_all_t_matrices(out: str) -> dict:
    """ 预扫描所有T矩阵文件建立快速索引 """
    mapping = defaultdict(dict)
    pattern = re.compile(r"Group_(\d+)_(\d+)_(\d+)_T_matrix\.csv$")
    
    for file in glob.glob(os.path.join(out, "Group_*_*_*_T_matrix.csv")):
        base = os.path.basename(file)
        match = pattern.match(base)
        if match:
            chrom, number = match.groups()[:2]
            mapping[(chrom, number)] = file
    return mapping

def process_group(args: tuple, t_matrix_cache: dict, t_matrix_mapping: dict, out: str) -> list:
    """ 处理单个组的核心函数 """
    (chrom, number), items = args
    results = []
    
    # 查找文件路径
    t_matrix_file = t_matrix_mapping.get((chrom, number), None)
    if not t_matrix_file:
        return [(idx, None, f"T matrix文件未找到，组: {chrom}_{number}") for idx, _, _, in items]

    # 读取数据 （带缓存）
    if t_matrix_file not in t_matrix_cache:
        try:
            t_matrix_cache[t_matrix_file] = pd.read_csv(t_matrix_file, index_col=0)
        except Exception as e:
            return [(idx, None, f"文件读取错误: {str(e)}") for idx, _, _ in items]
    
    df_t = t_matrix_cache[t_matrix_file]
    
    # 批量处理组内所有rSV
    for orig_idx, rSV_name, group_name in items:
        if rSV_name not in df_t.index:
            results.append((orig_idx, None, f"rSV数据缺失: {rSV_name}@{t_matrix_file}"))
            continue

        try:
            raw_gt = df_t.loc[rSV_name].values.tolist()
            converted_gt = list(map(convert_gt, raw_gt))  # 立即转换
            results.append((orig_idx, converted_gt, None))
        except Exception as e:
            results.append((orig_idx, None, f"GT处理异常: {str(e)}"))
    
    return results

def extract_vcf_sample(input_csv: str, output_gt: str, out: str, threads: int):
    """ 优化后的数据提取函数 """
    df_meta = pd.read_csv(input_csv)
    total = len(df_meta)
    
    # 初始化结果容器
    gt_matrix = [None] * total  # 保留原始顺序
    error_log = []

    # Stage 1: 预处理分组
    group_groups = defaultdict(list)
    for orig_idx, row in tqdm(df_meta.iterrows(), total=total, desc="Pre-processing groups", unit='rSV'):
        group_name = row["group_name"]
        match = re.fullmatch(r"Group_(\d+)_(\d+)_(\d+)_([rR][sS][vV]\d+)", group_name)
        
        if not match:
            error_log.append((orig_idx, f"Illegal group name: {group_name}"))
            continue

        chrom, number, _, rSV_name = match.groups()
        key = (chrom, number)
        group_groups[key].append((orig_idx, rSV_name, group_name))
    click.echo('Done')
    # Stage 2: 并行处理组
    t_matrix_cache = {}  # 文件路径到数据框的映射
    t_matrix_mapping = find_all_t_matrices(out)
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = []
        for group_key in group_groups:
            args = (group_key, group_groups[group_key])
            futures.append(executor.submit(
                process_group, args, t_matrix_cache, t_matrix_mapping, out
            ))
        group_count = len(futures)
        # 异步收集结果
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing group", unit='group'):
            try:
                group_results = future.result()
                for (idx, data, err) in group_results:
                    if err:
                        error_log.append((idx, err))
                    else:
                        gt_matrix[idx] = data
            except Exception as e:
                error_log.append((-1, f"Group process error: {str(e)}"))

    # Stage 3: 错误处理和输出
    valid_data = [data for data in gt_matrix if data is not None]
    gt_buffer = pd.DataFrame(valid_data)
    
    # 错误输出
    if error_log:
        error_df = pd.DataFrame(error_log, columns=["Original Line Number", "Error Details"])
        error_path = os.path.join(os.path.dirname(output_gt), "processing_errors.csv")
        error_df.to_csv(error_path, index=False)
        logger.warning(f"Warning: A total of {len(error_log)} errors were detected. See {error_path} for details.")

    logger.info(f"GT matrix successfully generated with {len(valid_data)} valid records.")
    return len(valid_data), group_count, gt_buffer
# 注：需要保留原convert_gt和其他辅助函数的实现

########生成rSV的VCF

def vcf_generate(sample_names: list, csv_file: str, gt_buffer: str, output_filled_vcf: str):
    # 读取 CSV 文件
    df = pd.read_csv(csv_file)
    # 将 GT DataFrame 转为 list of list（每一行是等价的）
    gt_list = gt_buffer.values.tolist()  # 每一行是个列表，等价于 split()
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
            logger.warning(f"Error processing row {index}: {e}")
            continue

    # 转换为 DataFrame
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    vcf_df = pd.DataFrame(vcf_data, columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"])

    # **按 CHROM 和 POS 进行排序**
    vcf_df.sort_values(by=["#CHROM", "POS"], inplace=True)

    """ 合并 VCF 数据和 GT 数据，直接写入最终文件 """
    with open(output_filled_vcf, "w") as out_vcf:
        # 写入头部信息
        out_vcf.write("##fileformat=VCFv4.2\n")
        out_vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        header = "\t".join(columns) + "\t" + "\t".join(sample_names)
        out_vcf.write(header + "\n")

        # 遍历输入数据的每一行
        for idx, vcf_row in enumerate(vcf_df.itertuples(index=False)):
            fields = list(map(str, vcf_row))          # VCF前9列
            gt_fields = list(map(str, gt_list[idx]))  # GT字段，转换为str
            fields[9:] = gt_fields                    # 替换掉后9列
            out_vcf.write("\t".join(fields) + "\n")


####检测错误GT
def detect_abnormal(out: str, config):
    # 获取当前目录下所有 T_matrix 文件
    t_matrix_files = [f for f in os.listdir(out) if re.match(r'Group_\d+_\d+_\d+_T_matrix\.csv', f)]
    # 统计信息列表
    stats_list = []
    # 添加前缀 'out/' 到每个文件名
    directory = out  # 'out' 是你要加前缀的目录
    t_matrix_files = [os.path.join(directory, file_name) for file_name in t_matrix_files]
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
    output_file = os.path.join(config.output_dir, "T_matrix_abnormal.csv")
    # 按照 abnormal_rate 降序排序
    df_stats_sorted = df_stats.sort_values(by="abnormal_rate", ascending=False)

    # 保存排序后的结果
    output_file = os.path.join(config.output_dir, "T_matrix_abnormal.csv")
    df_stats_sorted.to_csv(output_file, index=False)

    click.echo("Saved: T_matrix_abnormal.csv")

    # 计算所有 T_matrix 的汇总统计
    total_gt_0 = df_stats["number > 0"].sum()
    total_gt_2 = df_stats["number >2"].sum()
    total_abnormal_rate = total_gt_2 / total_gt_0 if total_gt_0 > 0 else 0

    # 生成 T_matrix_abnormal_all.csv
    df_total = pd.DataFrame([[total_gt_0, total_gt_2, total_abnormal_rate]], 
                            columns=["number > 0", "number >2", "abnormal_rate"])
    output_file = os.path.join(config.output_dir, "T_matrix_abnormal_all.csv")

    df_total.to_csv(output_file, index=False)
    click.echo("Saved: T_matrix_abnormal_all.csv")


