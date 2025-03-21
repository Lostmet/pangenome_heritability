import os
import glob
import re
import pandas as pd
import numpy as np
from tqdm import tqdm
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
from collections import defaultdict
from itertools import combinations
from ..utils.logging_utils import get_logger
logger = get_logger(__name__)


def parse_fasta_with_metadata(file_path: str):
    """
    Parse the given Fasta file, which is formatted like:
    
    >Group_2_2_50381
    C
    >Variant_2_2_50381_50381
    CTATGACTTTTGAATAAAAACCTAAAG
    >Variant_2_2_...
    ...
    >Group_2_3
    AA...
    >Variant_2_3_...

    Returns a data structure:
    {
      "Group_2_2_50381": {
          "reference": "C",
          "variants": [
              {
                  "chrom": 2,
                  "group": 2,
                  "start": 50381,
                  "end": 50381,
                  "sequence": "CTATGACTTTTGAATAAAAACCTAAAG"
              },
              ...
          ]
      },
      ...
    }
    """
    data = {}
    current_group = None
    current_reference = None

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">Group_"):
                # Parse group header, e.g., "Group_2_2"
                # Extract group identifier
                current_group = line[1:]  # Remove '>'
                data[current_group] = {
                    "reference": "",
                    "variants": []
                }
                current_reference = True  # The next line is the reference for this group
            elif line.startswith(">Variant_"):
                # If a Variant is encountered, parse its coordinate information
                variant_header = line[1:]  # e.g., "Variant_2_2_50381_50381"
                chrom, grp, start, end = parse_variant_header(variant_header)
                # Read the next line (or lines) as the sequence
                variant_seq = next(f).strip()  # Assume the variant always occupies one line

                if current_group is None:
                    raise ValueError("Found a Variant line before any Group line")

                data[current_group]["variants"].append({
                    "chrom": chrom,
                    "group": grp,
                    "start": start,
                    "end": end,
                    "sequence": variant_seq
                })
                current_reference = False
            else:
                # If it's neither a Group header nor a Variant header,
                # and current_reference=True, then this is the group's reference sequence
                if current_group and current_reference:
                    data[current_group]["reference"] += line
                else:
                    
                    pass
    #print(f"meta data: {data}")
    return data


def parse_variant_header(header: str):
    """
    Parse Variant_2_2_50381_50381
    Return (chrom, group, start, end)
    """
    # Assume the format is always correct: Variant_<chrom>_<group>_<start>_<end>
    pattern = r"^Variant_(\d+)_(\d+)_(\d+)_(\d+)$"
    match = re.match(pattern, header)
    if not match:
        raise ValueError(f"Variant header not in expected format: {header}")
    chrom = int(match.group(1))
    grp = int(match.group(2))
    start = int(match.group(3))
    end = int(match.group(4))
    return chrom, grp, start, end
def read_fasta_files(directory: str) -> Dict[str, List[Tuple[str, str]]]:
    """Read aligned FASTA files and handle missing/empty files."""
    pattern = os.path.join(directory, 'Group_*_*_aligned.fasta')
    file_paths = glob.glob(pattern)
    file_paths.sort()
    logger.info(f"Found {len(file_paths)} files matching the pattern.")
    
    fasta_contents = {}
    
    for file_path in file_paths:
        try:
            sequences = []
            with open(file_path, 'r') as file:
                current_seq = []
                seq_id = ''
                
                for line in file:
                    line = line.strip()
                    if line.startswith('>'):
                        if seq_id and current_seq:
                            sequences.append((seq_id, ''.join(current_seq)))
                        seq_id = line[1:]
                        current_seq = []
                    else:
                        current_seq.append(line)
                
                if seq_id and current_seq:
                    sequences.append((seq_id, ''.join(current_seq)))
            
            # blank file
            if not sequences:
                logger.warning(f"Empty aligned file detected: {file_path}")
                fasta_contents[os.path.basename(file_path)] = [('seq0', ''), ('seq1', '')]
            else:
                fasta_contents[os.path.basename(file_path)] = sequences
        
        except Exception as e:
            logger.error(f"Error processing file {file_path}: {str(e)}")
    
    return fasta_contents


def rSV_window_meta(sequence: str, start_offset: int) -> list:
    """
    Given a sequence and its start coordinate in the genome (start_offset),
    return [{pos: ..., rSV: ...}, ...], where pos is the absolute coordinate.
    """
    windows = []
    for i in range(len(sequence)):
        windows.append({
            'pos': start_offset + i,  # Calculate absolute coordinate
            'rSV': sequence[i:i + 1]
        })
    return windows


### <-- NEW OR MODIFIED CODE ###
def compare_rSVs_with_meta(ref_windows: List[Dict], var_windows: List[Dict]) -> Tuple[List[int], List[Dict]]:
    """
    Compare K-mer windows, treating '-' as a normal base for comparison.
    Deletion marks are only handled in the subsequent merging stage.
    
    Args:
        ref_windows: List of k-mer windows for the reference sequence
        var_windows: List of k-mer windows for the variant sequence
    
    Returns:
        Tuple[List[int], List[Dict]]: (difference array, metadata array)
    """
    if not ref_windows or not var_windows:
        logger.warning("Input windows are empty")
        return [0], [{'pos': 0, 'ref': '', 'alt': ''}]

    length = min(len(ref_windows), len(var_windows))
    
    diff_array = []
    meta_array = []

    for i in range(length): # 恢复对第一个的比对
        ref = ref_windows[i]
        var = var_windows[i]
        
        
        diff_val = 0 if ref["rSV"] == var["rSV"] else 1
        diff_array.append(diff_val)
        
        
        meta_array.append({
            'pos': ref["pos"],
            'ref': ref["rSV"],
            'alt': var["rSV"]
        })

    return diff_array, meta_array


def rSV_window(sequence: str, k: int = 4) -> List[str]:
    """
    Original simple version (if only a pure string list is needed).
    If rSV_window_meta is used, this function may no longer be needed.
    """
    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]


def compare_windows(ref_windows, var_windows):
    """
    Original compare_windows function.
    If compare_rSVs_with_meta is used, this function may no longer be needed.
    """
    if len(ref_windows) != len(var_windows):
        raise ValueError("Reference and variant windows must have the same length.")
    return [0 if ref == var else 1 for ref, var in zip(ref_windows, var_windows)]


def process_sequences(file_name: str, sequences: List[Tuple[str, str]], genome_metadata: dict ,has_insertion_dict: Dict) -> Dict:
    """
    Process sequence data, including simple seq0/seq1 cases and complex multi-sequence cases.
    """
    try:
       
        group_name = file_name.replace('_aligned.fasta', '').replace('_input.fasta', '')

        group_metadata = genome_metadata.get(group_name)
        #print(f"group_metadata: {group_metadata}")
        if not group_metadata:
            raise ValueError(f"Group {group_name} not found in genome metadata")
        
        results = []
        reference_seq = None
        reference_id = None            
    
        for seq_id, sequence in sequences:
            if seq_id.startswith('ref_') or seq_id == 'reference':
                reference_seq = sequence
                reference_id = seq_id
                break            
        
        if not reference_seq:
            reference_seq = sequences[0][1] 
            reference_id = sequences[0][0]            
        
        start_pos = group_metadata["variants"][0]["start"] if group_metadata["variants"] else 0
        
        # ref_windows：[{'pos': 906668, 'rSV': 'T'}, {'pos': 906669, 'rSV': '-'}]
        ref_windows = rSV_window_meta(reference_seq, start_pos)
        #print(f"ref_windows: {ref_windows}")
        number = 1
        for seq_id, sequence in sequences:
            if seq_id == reference_id:
                continue  
                
            var_windows = rSV_window_meta(sequence, start_pos)
            #print(f'var_windows: {type(var_windows)}')
            number += 1
            # 在这里进行逐个比对
            diff_array, meta_array = compare_rSVs_with_meta(ref_windows, var_windows)                
            results.append({
                'chromosome_group': group_name,
                'sequence_id': seq_id,
                'diff_array': diff_array,
                'meta_array': meta_array
            })

        return {
            'file_name': file_name,
            'results': results,
            'error': None
        }
        
    except Exception as e:
        logger.error(f"Error processing sequences in {file_name}: {str(e)}")
        return {'file_name': file_name, 'results': [], 'error': str(e)}




def process_fasta_files(
    directory: str,
    has_insertion_dict: dict, 
    genome_metadata: dict,  # Add genome_metadata parameter
    max_workers: Optional[int] = None,
    output_file: str = None,
    error_log: str = None
) -> Dict:
    """
    Process FASTA files and generate comparison results, integrating genome metadata.
    """
    try:
        if max_workers is None:
            max_workers = 10
        
        logger.info(f"Starting FASTA processing with {max_workers} workers")
        fasta_contents = read_fasta_files(directory)
        results = []
        errors = []
        
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(process_sequences, file_name, sequences, genome_metadata, has_insertion_dict): file_name
                for file_name, sequences in fasta_contents.items()
            }
            
            for future in tqdm(as_completed(futures), total=len(futures), desc="Generating rSVs"\
                               , unit='group'
                               ):
                result = future.result()
                if result['error']:
                    errors.append(f"Error in {result['file_name']}: {result['error']}")
                if result['results']:
                    results.extend(result['results'])
        
        if output_file:
            df = pd.DataFrame(results)
            df.to_csv(output_file, index=False)
            logger.info(f"Initial results saved to: {output_file}")
        
        if error_log and errors:
            with open(error_log, 'w') as f:
                for error in errors:
                    f.write(f"{error}\n")
            logger.warning(f"Errors were logged to: {error_log}")
        return {'processed': results, 'errors': errors}
    
    except Exception as e:
        logger.error(f"Error in process_fasta_files: {str(e)}")
        raise

def retain_changed_columns_group(rows: List[List[int]]) -> List[List[int]]:
    """
    rows: List of lists[int], each row is a diff_array: [0,1,0,1...].
    The function checks for changes across all rows in each column
    and returns a 2D list with only the columns where changes occurred.
    """
    if not rows:
        return []

    num_rows = len(rows)
    num_cols = len(rows[0])

    # Initialize the result list, keep the first column as is
    retained = [[row[0]] for row in rows]

    # Compare each column starting from the second one
    for col in range(1, num_cols):
        retain = False
        # Check each row to see if the current column differs from the previous column
        for row in range(num_rows):
            if row[col] != row[col - 1]:
                retain = True
                break

        # If there is a change, retain this column in all rows
        if retain:
            for row_i in range(num_rows):
                retained[row_i].append(rows[row_i][col])

    return retained

### <-- NEW OR MODIFIED CODE ###
def retain_changed_columns_group_with_index(rows: List[List[int]]) -> Tuple[List[List[int]], List[int]]:
    """
    Retain columns with changes and return their indices in the original rows[?].
    """
    if not rows:
        return [], []
    num_rows = len(rows)
    num_cols = len(rows[0])

    retained = [[] for _ in range(num_rows)]
    retained_indices = []

    # Retain the first column directly
    for r in range(num_rows):
        retained[r].append(rows[r][0])
    retained_indices.append(0)

    for col in range(1, num_cols):
        keep = False
        for row in rows:
            if row[col] != row[col - 1]:
                keep = True
                break
        if keep:
            for r in range(num_rows):
                retained[r].append(rows[r][col])
            retained_indices.append(col)

    return retained, retained_indices



def process_group(chrom, group, has_insertion_dict, cutoff):
    """
    处理每个chrom的组，返回处理结果。
    """
    '''
    # chrom:Group_2_1_1210539
    group:    
    chromosome_group        sequence_id    diff_array        meta_array
    0  Group_2_1_1210539        seq1  [0, 1, 1, 1, 1]  [{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'p...
    1  Group_2_1_1210539        seq2  [0, 1, 1, 1, 1]  [{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'p...
    2  Group_2_1_1210539        seq3  [0, 1, 1, 1, 1]  [{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'p...
    '''
    diff_arrays = []
    meta_arrays = []
    sequence_ids = []
    
    for _, row in group.iterrows():
        try:
            diff_array = eval(row['diff_array']) if isinstance(row['diff_array'], str) else row['diff_array']
            meta_array = eval(row['meta_array']) if isinstance(row['meta_array'], str) else row['meta_array']
            
            if not diff_array or not meta_array:
                logger.warning(f"Skipping empty arrays: chromosome_group={chrom}, sequence_id={row['sequence_id']}")
                continue
                
            diff_arrays.append(diff_array)
            meta_arrays.append(meta_array)
            sequence_ids.append(row['sequence_id'])
        except Exception as e:
            logger.warning(f"Error processing row data: {str(e)}, skipping this row. chromosome_group={chrom}, sequence_id={row['sequence_id']}")
            continue
    
    if len(diff_arrays) < 1 or len(meta_arrays) < 1:
        logger.warning(f"Group {chrom} does not have enough valid data for processing")
        return None
    
    # 使用 retain_changed_columns_group_with_meta 处理数据
    retained_diff, retained_meta, nrSV_set = retain_changed_columns_group_with_meta(diff_arrays, meta_arrays, chrom, has_insertion_dict, cutoff)
    if not retained_diff or not retained_meta or len(retained_diff) != len(sequence_ids):
        logger.warning(f"Invalid processing results for group {chrom}")
        return None
    
    #print(f'nrSV_set:{nrSV_set}')
    
    if nrSV_set == set():
        nrSV_merged_results = []
        nrSV_dict = {}
    else:
        nrSV_merged_results, nrSV_dict = nrSV_convert(chrom, nrSV_set)
    # 处理结果
    merged_results = []
    for i, seq_id in enumerate(sequence_ids):
        if i < len(retained_diff) and i < len(retained_meta):
            processed_meta = []
            for meta in retained_meta[i]:
                processed_meta.append({
                    'pos': meta['pos'],
                    'ref': meta['ref'],  
                    'alt': meta['alt']   
                })
            
            merged_results.append({
                'chromosome_group': chrom,
                'sequence_id': seq_id,
                'diff_array': retained_diff[i],
                'meta_array': processed_meta
            })
    #print('nrSV_merged_results: ', nrSV_merged_results)
    #print(f'merged_results: {merged_results}')
    return merged_results, nrSV_merged_results, nrSV_dict

import ast
from collections import defaultdict
from typing import List

def nrSV_convert(chrom: str, nrSV_set: set) -> List[dict]:
    merged_results = []
    nrSV_step = 0 # 储存一个步进变量
    nrSV_dict = {'group':chrom, 'SV':[]} # 储存一个dict用于后续vcf文件的sample便捷查找填补
    # 将字符串转换为字典
    sv_records = [ast.literal_eval(item) for item in nrSV_set]
    #print(f'sv_records{sv_records}')
    #[{'pos': 1210539, 'ref': 'G', 'alt': 'G', 'SV': 1}, {'pos': 1210542, 'ref': '-', 'alt': 'G', 'SV': 1}, {'pos': 1210541, 'ref': '-', 'alt': 'G', 'SV': 1}, {'pos': 1210543, 'ref': '-', 'alt': 'G', 'SV': 1}]
    # 根据 SV 值进行分组
    groups = defaultdict(list)
    for rec in sv_records:
        groups[rec['SV']].append(rec)
    
    # 对每个 SV 组，按 pos 排序后生成 meta_array，并生成一行结果
    for sv_value, records in sorted(groups.items()):
        #print(f'sv_value: {sv_value}')
        nrSV_dict['SV'].append(sv_value)
        # 按 pos 排序
        records.sort(key=lambda x: x['pos'])
        pos = records[0]['pos']
        ref = records[0]['ref']
        #print(f'sv_value:{sv_value}, records:{records}')
        #sv_value:1, records:[{'pos': 1210539, 'ref': 'G', 'alt': 'G', 'SV': 1}, {'pos': 1210541, 'ref': '-', 'alt': 'G', 'SV': 1}, {'pos': 1210542, 'ref': '-', 'alt': 'G', 'SV': 1}, {'pos': 1210543, 'ref': '-', 'alt': 'G', 'SV': 1}]
        merged = []  # 存放合并后的结果
        current = records[0].copy()   # 以第二个元素作为初始
        #print(f'records:{records}')
        if len(records) == 1:
            merged = records
        else:
            for rec in records[1:]:  # 从第2个开始遍历
                next_pos = current["pos"] + len(current["ref"])  # 计算合并后下一个期望位置
                
                if rec["pos"] == next_pos: 
                    # 合并插入突变
                    current["ref"] += rec["ref"]  # ref 继续扩展
                    current["alt"] += rec["alt"]  # alt 也累加
                else:
                    # 不能合并，存储当前合并块，并重新开始
                    merged.append(current)
                    current = rec.copy()

            # 别忘了最后一个片段
            merged.append(current)
        #print(f'merged:{merged}')
        # 合并所有记录的 meta 信息
        meta_array = [{'pos': pos, 'ref': ref+rec['ref'].replace('-',''), 'alt': ref+rec['alt'].replace('-','')} for rec in merged[1:]]
        nrSV_count = len(meta_array)
        for i in range(nrSV_count):
            merged_results.append({
                'chromosome_group': chrom + f'_nrSV{i+1+nrSV_step}',
                'sequence_id': f"seq{sv_value}",  # SV 值对应序列号，如 SV==1 -> seq1
                'meta_array': meta_array[i]
            })
        nrSV_step += nrSV_count
    
    # 如果有多个组，按序号排序（这里 sequence_id 格式为 "seq{number}"）
    merged_results.sort(key=lambda x: int(x['chromosome_group'].split('_')[4][4:]))
    # print(f'merged_results: {merged_results}')
    # [{group:, seqID:, meta_array:[{},{},...,{}]}]
    return merged_results, nrSV_dict

import pysam
def nrSV_vcf_generate(nrSV_csv: str, config, nrSV_list, nSV_name: str):
    df = pd.read_csv(nrSV_csv) 
    nrSV_count = df.shape[0]
    output_cache = {}  # 存储结果的字典
    # 拿到output的vcf的名字
    output_vcf = os.path.join(config.output_dir, nSV_name)
    # 读取原本vcf的信息
    with pysam.VariantFile(config.vcf_file, "r") as vcf_in:
        for item in nrSV_list:
            parts = item['group'].split('_')
    
            chrom = parts[1]
            pos = int(parts[3])
            sv_indexes = item['SV']
            # 初始化存储该 group 的基因型数据
            output_cache[item['group']] = {}
        
            # 遍历 VCF 该位点的所有变异
            sv_records = []
            for record in vcf_in.fetch(chrom, pos - 1, pos):  # pos-1 因为 VCF 使用 1-based 但 fetch 使用 0-based
                if record.pos == pos:
                    sv_records.append(record)
        
            # 提取 SV 索引对应的变异
            for sv_idx in sv_indexes:
                if sv_idx <= len(sv_records):  # 确保索引有效
                    sv_record = sv_records[sv_idx - 1]  # SV 索引是从 1 开始的，而 Python list 是 0-based

                    # 获取所有 sample 的 GT 信息
                    gt_data = {sample: sv_record.samples[sample]['GT'] for sample in sv_record.samples}

                    # 存储基因型数据
                    output_cache[item['group']][f"SV{sv_idx}"] = gt_data 
                    # print(f'idx:{sv_idx},gt:{gt_data}')
    # 获取 output_cache 里第一个可用的 group_name 和 SV_x 作为 sample_names 的来源
    sample_names = []

    if output_cache:  # 确保缓存不为空
        first_group = next(iter(output_cache))  # 获取第一个 group_name
        first_sv = next(iter(output_cache[first_group])) if output_cache[first_group] else None  # 获取第一个 SV_x

        if first_sv:
            sample_names = list(output_cache[first_group][first_sv].keys())  # 获取 sample 名称列表
    # 写入输出
    with open(output_vcf, "a") as f_out:
        for _, row in df.iterrows():
            group_name_parts = row['chromosome_group'].split('_')
            group_name = '_'.join(group_name_parts[:-1])
            nrSV_id = row['chromosome_group']
            chrom = group_name_parts[1]
            meta = ast.literal_eval(row['meta_array'])
            sv_id = int(row['sequence_id'][3:])
            ref = meta['ref']
            alt = meta['alt']
            pos = meta['pos']
            try:
                gt_info = output_cache[group_name][f'SV{sv_id}']
            except KeyError as e:
                print(f"KeyError: {e} not found in output_cache")
                gt_info = {}
            gt_values = ["/".join(map(str, gt_info.get(sample, ('.', '.')))) for sample in sample_names]
            f_out.write(f"{chrom}\t{pos}\t{nrSV_id}\t{ref}\t{alt}\t.\t.\ttype=nrSV\tGT\t" + "\t".join(gt_values) + "\n")
    return nrSV_count

def process_and_merge_results(results, output_csv: str, threads: str, \
                              has_insertion_dict: dict, cutoff: float, nrSV_csv: str):
    """
    Process and merge k-mer window results, removing '-' characters in the final results.
    """
    # has_insertion_dict:{'Group_2_1_1210539': True}
    try:
        df = pd.DataFrame(results['processed'])
    except Exception as e:
        logger.error(f"Error saving results to CSV: {str(e)}")
        raise
        
    try:
        if df.empty:
            logger.warning("Input file is empty")
            return
        
        # 根据chromosome_group分组
        grouped = df.groupby('chromosome_group')
        
        # 计算进度条总任务数：所有组的数量
        total_groups = len(grouped)
        
        merged_results = []
        nrSV_merged_results = []
        nrSV_list = []
        # 使用进度条
        with tqdm(total=total_groups, desc="c. Merging rSVs", unit="group") as pbar:
            # 使用ProcessPoolExecutor来并行处理每个组
            with ProcessPoolExecutor(max_workers=threads) as executor:
                future_to_group = {executor.submit(process_group, chrom, group, has_insertion_dict, cutoff): chrom for chrom, group in grouped}
                # 遍历每个任务的结果
            for future in as_completed(future_to_group):
                try:
                    result, nrSV_result, nrSV_dict = future.result()
                    if result:
                        merged_results.extend(result)
                    if nrSV_result:
                        nrSV_merged_results.extend(nrSV_result)
                    if nrSV_dict:
                        nrSV_list.append(nrSV_dict)
                except Exception as e:
                    group = future_to_group[future]
                    print(f"❌ Error in processing group {group}: {e}")
                    print(f"💡 Tip: Please check whether the sequences in group {group} are normalized properly.")
                pbar.update(1)

        if not merged_results:
            logger.warning("No valid merged results")
            return
        
        # 创建新的DataFrame并保存
        merged_df = pd.DataFrame(merged_results)
        merged_df.to_csv(output_csv, index=False)
        logger.info(f"Merged results of rSV saved to: {output_csv}")
        # 生成nrSV_results.csv
        nrSV_merged_df = pd.DataFrame(nrSV_merged_results)
        nrSV_merged_df.to_csv(nrSV_csv, index=False)
        logger.info(f"Merged results of nrSV saved to: {nrSV_csv}")
    except Exception as e:
        logger.error(f"Error during processing: {str(e)}")
        raise
    return nrSV_list
    # print(f'merge输入threads:{threads}')

def process_chromosome_groups(input_csv: str, output_csv: str) -> None:
    """
    Process data for each chromosome_group, retaining arrays for all positions.
    """
    data = defaultdict(list)
    
    # Read CSV file
    with open(input_csv, newline='', encoding='utf-8') as csvfile:
        reader = pd.read_csv(csvfile)
        
        for _, row in reader.iterrows():
            chromosome_group = row['chromosome_group']
            diff_array = eval(row['diff_array']) if isinstance(row['diff_array'], str) else row['diff_array']
            meta_array = eval(row['meta_array']) if isinstance(row['meta_array'], str) else row['meta_array']
            
            # Directly add original data to processed_data
            #processed_data.append({
            #    'chromosome_group': chromosome_group,
            #    'sequence_id': row['sequence_id'],
            #    'diff_array': diff_array,  # Use original diff_array
            #    'meta_array': meta_array   # Use original meta_array
            #})
    
    # Write processed CSV
    #processed_df = pd.DataFrame(processed_data)
    # Convert arrays to strings for storage
    #processed_df['diff_array'] = processed_df['diff_array'].apply(str)
    #processed_df['meta_array'] = processed_df['meta_array'].apply(str)
    #processed_df.to_csv(output_csv, index=False, encoding='utf-8')
    logger.info(f"Processed results saved to: {output_csv}")


def save_rSV_results_to_csv(results: Dict, output_file: str) -> None:
    """Save rSV comparison results to a CSV file."""
    try:
        if results['processed']:
            df = pd.DataFrame(results['processed'])
            df.to_csv(output_file, index=False)
            logger.info(f"Results saved to: {output_file}")
        else:
            logger.warning("No results to save")
    except Exception as e:
        logger.error(f"Error saving results to CSV: {str(e)}")
        raise


### <-- NEW OR MODIFIED CODE ###
def explode_final_results(input_csv: str, output_csv: str) -> None:
    """
    To expand meta_array into multiple rows, one k-mer per row, and write pos/ref/alt/diff information
    to the table, use this function as an example.
    """
    df = pd.read_csv(input_csv)
    # Assuming there are 'chromosome_group', 'sequence_id', 'diff_array', 'meta_array'
    all_rows = []
    for _, row in df.iterrows():
        cg = row['chromosome_group']
        seq_id = row['sequence_id']
        diff_array = eval(row['diff_array']) if isinstance(row['diff_array'], str) else row['diff_array']
        meta_array = eval(row['meta_array']) if isinstance(row['meta_array'], str) else row['meta_array']

        for diff_val, meta in zip(diff_array, meta_array):
            # meta: {'pos':..., 'ref':..., 'alt':...}
            out_row = {
                'chromosome_group': cg,
                'sequence_id': seq_id,
                'pos': meta['pos'],
                'ref': meta['ref'],
                'alt': meta['alt'],
                'difference': diff_val
            }
            all_rows.append(out_row)

    df_out = pd.DataFrame(all_rows)
    df_out.to_csv(output_csv, index=False)
    logger.info(f"Exploded results saved to: {output_csv}")

def retain_changed_columns_group_with_meta(
    rows: List[List[int]], 
    meta_rows: List[List[Dict]],
    chrom,
    has_insertion_dict,
    cutoff # 注意是个float，是百分比
) -> Tuple[List[List[int]], List[List[Dict]]]: # 合并函数
    """
    Process difference arrays and metadata arrays, merging identical columns starting from the first column.
    
    Args:
        rows: List of difference arrays, each row represents the alignment result of a sequence
        meta_rows: List of metadata, corresponding to the difference arrays
    
    Returns:
        Tuple[List[List[int]], List[List[Dict]]]: (Processed difference arrays, Processed metadata arrays)
    """
    if not rows or not meta_rows:
        return [], []
    
    if len(rows) != len(meta_rows):
        raise ValueError(f"Mismatch between the lengths of difference arrays and metadata arrays, meta_data_pos:{meta_rows[0][0]['pos']}")
    # 调试print meta_rows
    # print(f"Mismatch between the lengths of difference arrays and metadata arrays, meta_data_pos:{meta_rows[0][0]['pos']}")
    # print(f'正在处理:{chrom}')
    num_rows = len(rows)
    num_cols = len(rows[0])
    #print(f'num_rows:{num_rows}, num_cols:{num_cols}')

    # Initialize result lists, 这个初始化看上去是把第一个东西先扔上去
    retained_diff = [[rows[i][0]] for i in range(num_rows)]
    retained_meta = [[meta_rows[i][0]] for i in range(num_rows)]
    #print(f'diff:{rows}, meta:{meta_rows}\nretained diff:{retained_diff}, meta: {retained_meta}')
    # diff:[[0, 1, 1, 1, 1, 1, 1], [1, 1, 0, 0, 0, 0, 0]]
    # retained diff:[[0], [1]]，保留了每组第一个
    current_col = 0  # Current column index being processed
    
    # Compare starting from the second column
    for col in range(1, num_cols):
        is_same = True
        # Check if all rows in the current column are the same as the previous column
        for row_i in range(num_rows):
            if rows[row_i][col] != rows[row_i][col-1]:
                is_same = False
                break
        
        if is_same:
            # If the current column is the same as the previous column, merge metadata
            for row_i in range(num_rows):
                curr_meta = meta_rows[row_i][col]
                prev_meta = retained_meta[row_i][-1]
                
                # Merge metadata, adding new bases
                merged_meta = {
                    'pos': prev_meta['pos'],
                    'ref': prev_meta['ref'] + curr_meta['ref'][-1],
                    'alt': prev_meta['alt'] + curr_meta['alt'][-1]
                }
                retained_meta[row_i][-1] = merged_meta
        else:
            # If different, add new column
            for row_i in range(num_rows):
                retained_diff[row_i].append(rows[row_i][col])
                retained_meta[row_i].append(meta_rows[row_i][col])
            current_col += 1

    nrSV_set = set()
    new_rSV_count = 0 # 记录多少新的rSV由于重叠但不相同而拆开
    ## 下面开始写对poly_ins中可能的poly_alt的情况的diff和meta的改动代码
    # 首先判断一下是不是poly_ins组，不是就直接跳过了
    if has_insertion_dict[chrom]:# 也就是if true，{'Group_2_1_1210539': True} 
        sv_count = len(retained_diff) # 先确定循环的次数，也就是sv的个数
        rSV_count = len(retained_diff[0]) # 确定初始rSV的个数
        rSV_index_small = 0 # 偏移值
        poly_pos = 0 # 用于储存pos的index的步进补偿（用于后续的标准化）

        for i in range(1, rSV_count): # 遍历现在的rSV（第一个凑数的不用看所以从1开始）

            rSV_index = i + rSV_index_small #####用于补偿由于插入带来的i偏移
           
            #print(f'i:{i},rSV_index:{rSV_index}') 
            alt_list = [] # 创建一个列表，用于无差别存储alt以获取多alt的重复情况
            alt_set = set() # 创建一个集合，来看是否是多alt
           
            for j in range(sv_count): # 对固定rSV（列）的每行进行遍历
                if retained_diff[j][rSV_index] == 1 \
                    and retained_meta[j][rSV_index]['alt'][0] != "-" : # 如果对应的第i个rSV的SV（遍历j）为1，且alt不为空，则添加到集合里
                    #start_count = len(alt_set) 
                    alt_origin = retained_meta[j][rSV_index]['alt']
                    alt_set.add(alt_origin) #  'alt': 'GGGG'
                    alt_list.append(alt_origin) # 存储一个集合和一个列表
                    #end_count = len(alt_set)
                    #if start_count != end_count:
                        #alt_index_list.append(j) # 获取对应的alt的索引
            
            #print(f'alt_index_list:{alt_index_list}')
            #print(f'alt_set:{alt_set}')
            poly_alt_count = len(alt_set) # 查一下有多少个alt
                
            if poly_alt_count <= 1: # 那就说明没有多alt的情况了，只需要处理pos
                for j in range(sv_count): # 对对应rSV的所有SV的pos进行处理
                    retained_meta[j][rSV_index]['pos'] += poly_pos
            # 在这里试着添加alt相似阈值，逻辑：俩俩比较，只要有一个满足阈值则保留，一个都不满足则剔除
            

            else:  # 不然需要对diff_array和meta_array进行填充处理
                '''首先对重叠的alt进行提取，也就是alt_set:{'G', 'A'}类的通过filter之后我们不把这个当成新的rSV'''
                ## 下面添加阈值
                alt_set_filtered = filter_similar_sequences(alt_list, cutoff) # 输入一个阈值
                #print(f'alt_list:{alt_list} alt_set_filtered:{alt_set_filtered}')
                ###########################################
                # step 0. 如果一个没保留，直接把列删了。
                ###########################################
                if alt_set_filtered == set(): # 说明一个都没保留, 需要删除该rSV的diff_array和meta_array
                    
                    poly_pos -= len(retained_meta[0][rSV_index]['ref'])
                    for j in range(sv_count):
                        retained_diff[j].pop(rSV_index) # 删除diff
                        nrSV_meta = retained_meta[j][rSV_index] # {ref:xxx, alt:xxx, pos:xxx}
                        nrSV_meta['SV'] = j+1 # {ref:xxx, alt:xxx, pos:xxx, SV:xxx}
                        nrSV_set.add(str(nrSV_meta)) # "{ref:xxx, alt:xxx, pos:xxx, SV:xxx}"
                        meta0 = retained_meta[j][0] # {ref: A, alt : A}, 最前面用于定位的
                        meta0['SV'] = j+1
                        nrSV_set.add(str(meta0))


                        retained_meta[j].pop(rSV_index) # 删除meta
                    # 去步进
                    rSV_index_small -= 1 

                else:
                    
                    new_rSV_count += len(alt_set_filtered) - 1
                    #print(f'alt_set_filtered:{alt_set_filtered}')
                    # 得到的set：从 ["AAAAAAA", "AAAAAAA", "AAAAAAG", "GGGGGGG", "GGCACAG"]
                    # 阈值为90%得到 {'AAAAAAA'}
                    #############################################
                    # step 1: 先拿到了保留重叠项的alt_set_filtered
                    #############################################
                    #print('开始处理多alt')
                    poly_alt = 0 # 用于储存alt的index的步进
                    alt_index_list = [] # 创建一个列表，得到后续哪个SV是poly_alt的，比如[2,3]对应了alt_set的的{'AAA','CCC'}
                    meta_list = [] # 储存一个不被index污染的metalist
                    alt_nrSV_index_set = set() # {1,2...} nrSV的对应的j
                    #############################################
                    # step 2: 处理非空未重叠项（nrSV）
                    #############################################
                    for j in range(sv_count): ###先遍历SV得到被排除的nrSV的meta
                        alt = retained_meta[j][rSV_index]['alt'] # 比如说GGG
                        meta = retained_meta[j][rSV_index] # 缓存一下meta                   
                        if alt not in alt_set_filtered and alt[0] != "-":# 如果alt对不上，并且非空，那么说明没有重复也没达到重叠的阈值
                            meta['SV'] = j+1 # 给meta对应的字典多加一个SV的对应编号（或许不需要+1，再看）
                            # 现在meta是{'pos':xxxx, 'ref':xxxxx, 'alt':xxxxx, 'SV':x}
                            nrSV_set.add(str(meta)) # 嵌套，但是天然去重
                            meta0 = retained_meta[j][0] # {ref: A, alt : A}, 最前面用于定位的
                            meta0['SV'] = j+1  
                            nrSV_set.add(str(meta0))                          
                            # set：{{ref:xxx, alt:xxx .....}}
                            # 接着把原始位置的diff搞成0, 表示这个无关了
                            retained_diff[j][rSV_index] = 0
                            alt_nrSV_index_set.add(j)
                            #print(f'meta_nrSV_set:{meta_nrSV_set}')
                        #nrSV_meta_list_all.append(meta_nrSV_set)

                    for alt_filtered in alt_set_filtered:  # 遍历每个候选alt，比如说 AAA
                        found = False  # 用于标记是否找到对应的alt
                        for j in range(sv_count):  # 遍历所有 SV
                            #下面试着添加一下，nrSV的内容
                            if retained_meta[j][rSV_index]['alt'] == alt_filtered and not found:  # 第一次找到对应 alt
                                alt_index_list.append(j)  # 记录 sv 位置
                                found = True  # 标记找到，保证只会找到的第一次来到这个位置
                                meta = retained_meta[j][rSV_index] # 缓存一下meta
                                break

                        meta['pos'] += poly_pos
                        #last_length = len(meta['alt']) # 储存一个last_length用于步进削减
                        #poly_pos += last_length # 储存步进
                        #print(f"meta:{meta},len:{len(meta['alt'])}")
                        meta_list.append(meta)        
                    #print(f'meta_list:{meta_list},alt_index_list:{alt_index_list}')
                    _count_dict = {} # 初始化占位符count的字典用于查验
                    
                    for index, alt_index in enumerate(alt_index_list): # 遍历alt的index，比如对第一个SV进行rSV的扩张变换
                        meta = meta_list[index] # 这里应该是一个meta={}，表示对应的alt的meta
                        alt_to_compare = retained_meta[alt_index][rSV_index+ poly_alt]['alt'] # 列表中进行比较的
                        #print(f'对第alt_index={alt_index}的循环，meta为{meta}')                 
                        for j in range(sv_count): # 对某一个具体的变化进行逐个遍历
                            # 记录alt出现占位符的次数
                            #print(f'there, j={j}, index={rSV_index+poly_alt},retained_meta[j]={retained_meta[j]},')
                            if rSV_index + poly_alt < len(retained_meta[j]) :
                                alt = retained_meta[j][rSV_index+ poly_alt]['alt'] # 
                            else:
                                alt = None
                            #print(f'here, j={j}, alt=:{alt}')
                            alt_del = retained_meta[j][rSV_index+ poly_alt-index]['alt'] #对del进行去步进                             
                            #print(f"{alt}对比{alt_to_compare}")
                            if alt_del[0] == "-" or j in alt_nrSV_index_set: # 如果是占位符或者因为无足够重叠要被删掉的nrSV，特殊处理
                                if f'SV{j}' not in _count_dict: # 如果该SV还没出现过
                                    _count_dict[f'SV{j}'] = j
                                    #print(f'1_count_dict:{_count_dict}')
                                    retained_meta[j][rSV_index + poly_alt]['pos'] = meta['pos'] # 只改pos
                                elif _count_dict[f'SV{j}'] == j: # 该SV出现过一次，则加meta
                                    #print(f'2_count_dict:{_count_dict}')
                                    retained_diff[j][rSV_index + poly_alt] =  0 # 就对对应rSV位置的前面加一个0(insert是索引前的位置插入)
                                    retained_meta[j][rSV_index + poly_alt] = meta  
                                
                            elif alt != alt_to_compare: # 如果并非该alt变异，就对其diff和meta进行扩充

                                #print(f'到这儿3, j={j}')                    
                                retained_diff[j][rSV_index + poly_alt] =  0 # 就对对应rSV位置的前面加一个0(insert是索引前的位置插入)
                                retained_meta[j][rSV_index + poly_alt] = meta
                                #print(f'第{j}个SV，在{rSV_index+poly_alt-1}的位置插入了{meta}')
                            
                            else:                               
                                retained_meta[j][rSV_index + poly_alt] = meta                                                     
                                #print(f'我没通过,meta是{meta},j={j},i={i},rSV_index={rSV_index}')
                        #poly_alt += 1 # 储存步进
                    
                    #rSV_index_small += len(meta_list)-1
                    #poly_pos -=  last_length# 还得减去poly的偏移值

    return retained_diff, retained_meta, nrSV_set 
# 输出的是一个组的比如：
##retained_diff: meta：
# [[0, 1], [[{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '----', 'alt': 'GGGG'}],
# [0, 1],  [{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '----', 'alt': 'AAAA'}], 
# [0, 1]], [{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '----', 'alt': 'ATAA'}]]

##retained_diff: meta：
# [[0, 1],    [[{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '----', 'alt': 'GGGG'}],
# [0, 0, 1],  [{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '----', 'alt': 'GGGG'}, {'pos': 1210540, 'ref': '----', 'alt': 'AAAA'}], 
# [0, 0, 1]], [{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '----', 'alt': 'GGGG'}, {'pos': 1210540, 'ref': '----', 'alt': 'ATAA'}]]

'''
i = 1
diff:
[[0, 0, 1, 1, 1], 
[0, 0, 1, 1, 0], 
[0, 1, 1, 0, 0]]
meta:
[[{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '-', 'alt': '-'}, {'pos': 1210541, 'ref': '--', 'alt': 'GG'}, {'pos': 1210543, 'ref': '-', 'alt': 'G'}, {'pos': 1210544, 'ref': '-', 'alt': 'G'}], 
[{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '-', 'alt': '-'}, {'pos': 1210541, 'ref': '--', 'alt': 'AA'}, {'pos': 1210543, 'ref': '-', 'alt': 'A'}, {'pos': 1210544, 'ref': '-', 'alt': '-'}], 
[{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '-', 'alt': 'T'}, {'pos': 1210541, 'ref': '--', 'alt': 'AA'}, {'pos': 1210543, 'ref': '-', 'alt': '-'}, {'pos': 1210544, 'ref': '-', 'alt': '-'}]]

i = 2
diff:
[[0, 0, 1, 0, 1, 1], 
[0, 0, 0, 1, 1, 0], 
[0, 1, 0, 0, 1, 0, 0]]
meta:
[[{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '-', 'alt': '-'}, {'pos': 1210541, 'ref': '--', 'alt': 'GG'}, {'pos': 1210543, 'ref': '--', 'alt': 'AA'}, {'pos': 1210543, 'ref': '-', 'alt': 'G'}, {'pos': 1210544, 'ref': '-', 'alt': 'G'}], 
[{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '-', 'alt': '-'}, {'pos': 1210541, 'ref': '--', 'alt': 'GG'}, {'pos': 1210543, 'ref': '--', 'alt': 'AA'}, {'pos': 1210543, 'ref': '-', 'alt': 'A'}, {'pos': 1210544, 'ref': '-', 'alt': '-'}], 
[{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '-', 'alt': 'T'}, {'pos': 1210541, 'ref': '--', 'alt': 'GG'}, {'pos': 1210543, 'ref': '--', 'alt': 'AA'}, {'pos': 1210541, 'ref': '--', 'alt': 'AA'}, {'pos': 1210543, 'ref': '-', 'alt': '-'}, {'pos': 1210544, 'ref': '-', 'alt': '-'}]]

i = 2改
diff:
[[0, 0, 1, 1, 1], 
[0, 0, 0, 1, 1, 0], 
[0, 1, 0, 1, 0, 0]]
mata:
[[{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '-', 'alt': '-'}, {'pos': 1210541, 'ref': '--', 'alt': 'GG'}, {'pos': 1210543, 'ref': '-', 'alt': 'G'}, {'pos': 1210544, 'ref': '-', 'alt': 'G'}], 
[{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '-', 'alt': '-'}, {'pos': 1210541, 'ref': '--', 'alt': 'GG'}, {'pos': 1210543, 'ref': '--', 'alt': 'AA'}, {'pos': 1210543, 'ref': '-', 'alt': 'A'}, {'pos': 1210544, 'ref': '-', 'alt': '-'}], 
[{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '-', 'alt': 'T'}, {'pos': 1210541, 'ref': '--', 'alt': 'GG'}, {'pos': 1210541, 'ref': '--', 'alt': 'AA'}, {'pos': 1210543, 'ref': '-', 'alt': '-'}, {'pos': 1210544, 'ref': '-', 'alt': '-'}]]

i =2 改i的步进，对了
diff:
[[0, 0, 1, 0, 1, 1], 
[0, 0, 0, 1, 1, 0], 
[0, 1, 0, 1, 0, 0]]
[[{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '-', 'alt': '-'}, {'pos': 1210541, 'ref': '--', 'alt': 'GG'}, {'pos': 1210543, 'ref': '--', 'alt': 'AA'}, {'pos': 1210543, 'ref': '-', 'alt': 'G'}, {'pos': 1210544, 'ref': '-', 'alt': 'G'}], 
[{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '-', 'alt': '-'}, {'pos': 1210541, 'ref': '--', 'alt': 'GG'}, {'pos': 1210543, 'ref': '--', 'alt': 'AA'}, {'pos': 1210543, 'ref': '-', 'alt': 'A'}, {'pos': 1210544, 'ref': '-', 'alt': '-'}], 
[{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '-', 'alt': 'T'}, {'pos': 1210541, 'ref': '--', 'alt': 'GG'}, {'pos': 1210541, 'ref': '--', 'alt': 'AA'}, {'pos': 1210543, 'ref': '-', 'alt': '-'}, {'pos': 1210544, 'ref': '-', 'alt': '-'}]]

i =3
diff:
[[0, 0, 1, 0, 1, 0, 1], 
[0, 0, 0, 1, 0, 1, 0], 
[0, 1, 0, 1, 0, 0, 0, 0]]
mata:
[[{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '-', 'alt': '-'}, {'pos': 1210541, 'ref': '--', 'alt': 'GG'}, {'pos': 1210543, 'ref': '--', 'alt': 'AA'}, {'pos': 1210543, 'ref': '-', 'alt': 'G'}, {'pos': 1210544, 'ref': '-', 'alt': 'A'}, {'pos': 1210544, 'ref': '-', 'alt': 'G'}], 
[{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '-', 'alt': '-'}, {'pos': 1210541, 'ref': '--', 'alt': 'GG'}, {'pos': 1210543, 'ref': '--', 'alt': 'AA'}, {'pos': 1210543, 'ref': '-', 'alt': 'G'}, {'pos': 1210544, 'ref': '-', 'alt': 'A'}, {'pos': 1210544, 'ref': '-', 'alt': '-'}], 
[{'pos': 1210539, 'ref': 'G', 'alt': 'G'}, {'pos': 1210540, 'ref': '-', 'alt': 'T'}, {'pos': 1210541, 'ref': '--', 'alt': 'GG'}, {'pos': 1210541, 'ref': '--', 'alt': 'AA'}, {'pos': 1210543, 'ref': '-', 'alt': 'G'}, {'pos': 1210544, 'ref': '-', 'alt': 'A'}, {'pos': 1210543, 'ref': '-', 'alt': '-'}, {'pos': 1210544, 'ref': '-', 'alt': '-'}]]
'''
 # 输入dna_set和threshold，输出新的set
def hamming_similarity(seq1, seq2):
    """计算两个等长 DNA 序列的 Hamming 相似度"""
    assert len(seq1) == len(seq2), "序列长度必须相同"
    return sum(a == b for a, b in zip(seq1, seq2)) / len(seq1)

# 现在所有重复的DNA序列也会被保留下来了
def filter_similar_sequences(dna_list, threshold):
    """保留相似度大于等于 threshold 的 DNA 序列"""
    #dna_list = list(dna_set)  # 转为列表
    n = len(dna_list)
    similarity_matrix = np.zeros((n, n))

    # 计算相似度矩阵
    for i, j in combinations(range(n), 2):
        similarity_matrix[i, j] = similarity_matrix[j, i] = hamming_similarity(dna_list[i], dna_list[j])

    # 统计每个序列与其他序列的最大相似度
    max_similarities = similarity_matrix.max(axis=1)

    # 只保留至少有一个相似度 ≥ threshold 的序列
    filtered_sequences = {dna_list[i] for i in range(n) if max_similarities[i] >= threshold}

    return filtered_sequences
