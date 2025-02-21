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


def kmer_window_meta(sequence: str, start_offset: int) -> list:
    """
    Given a sequence and its start coordinate in the genome (start_offset),
    return [{pos: ..., kmer: ...}, ...], where pos is the absolute coordinate.
    """
    windows = []
    for i in range(len(sequence)):
        windows.append({
            'pos': start_offset + i,  # Calculate absolute coordinate
            'kmer': sequence[i:i + 1]
        })
    return windows


### <-- NEW OR MODIFIED CODE ###
def compare_kmers_with_meta(ref_windows: List[Dict], var_windows: List[Dict]) -> Tuple[List[int], List[Dict]]:
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

    for i in range(length - 1): # 去掉对第一个的比对
        ref = ref_windows[i + 1]
        var = var_windows[i + 1]
        
        
        diff_val = 0 if ref["kmer"] == var["kmer"] else 1
        diff_array.append(diff_val)
        
        
        meta_array.append({
            'pos': ref["pos"],
            'ref': ref["kmer"],
            'alt': var["kmer"]
        })

    return diff_array, meta_array


def kmer_window(sequence: str, k: int = 4) -> List[str]:
    """
    Original simple version (if only a pure string list is needed).
    If kmer_window_meta is used, this function may no longer be needed.
    """
    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]


def compare_windows(ref_windows, var_windows):
    """
    Original compare_windows function.
    If compare_kmers_with_meta is used, this function may no longer be needed.
    """
    if len(ref_windows) != len(var_windows):
        raise ValueError("Reference and variant windows must have the same length.")
    return [0 if ref == var else 1 for ref, var in zip(ref_windows, var_windows)]


def process_sequences(file_name: str, sequences: List[Tuple[str, str]], genome_metadata: dict ,has_insertion_dict: Dict) -> Dict:
    """
    Process sequence data, including simple seq0/seq1 cases and complex multi-sequence cases.
    
    Args:
        file_name: The name of the current file
        sequences: List of sequence data [(seq_id, sequence), ...]
        genome_metadata: Genome metadata
        k: k-mer window size
    
    Returns:
        Dict: A dictionary containing the comparison results
    """
    #print(f"file_name: {file_name}, sequences:{sequences}")
    try:
       
        group_name = file_name.replace('_aligned.fasta', '').replace('_input.fasta', '')

        if not has_insertion_dict[group_name]:
            #print(f'成功通过该检查点，也就是not false, group为{group_name}')

            # 对未重叠组的操作，似乎是冗余的，未重叠组只需要保留就可以了，不过多一个检索也不是问题，但是没必要接着做
            if len(sequences) == 2 and sequences[0][0] == 'seq0' and sequences[1][0] == 'seq1':
                ref_seq = sequences[0][1]  # seq0
                var_seq = sequences[1][1]  # seq1
            
                group_metadata = genome_metadata.get(group_name)
                if not group_metadata or not group_metadata["variants"]:
                    logger.warning(f"No metadata found for group {group_name}")
                    start_pos = 0
                else:
                    start_pos = group_metadata["variants"][0]["start"]
                
                diff_array = []
                meta_array = []                
                
                for i in range(len(ref_seq)):
                    ref_kmer = ref_seq[i:i+1]
                    var_kmer = var_seq[i:i+1]                    
                    
                    diff_array.append(1)                   
                    
                    meta_array.append({
                        'pos': start_pos + i,
                        'ref': ref_kmer,  
                        'alt': var_kmer  
                    })
                
                return {
                    'file_name': file_name,
                    'results': [{
                        'chromosome_group': group_name,
                        'sequence_id': 'seq1',
                        'diff_array': diff_array,
                        'meta_array': meta_array
                    }],
                    'error': None
                }
            
        
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
            
            # ref_windows：[{'pos': 906668, 'kmer': 'T'}, {'pos': 906669, 'kmer': '-'}]
            ref_windows = kmer_window_meta(reference_seq, start_pos)
            #print(f"ref_windows: {ref_windows}")
            number = 1
            for seq_id, sequence in sequences:
                if seq_id == reference_id:
                    continue  
                    
                var_windows = kmer_window_meta(sequence, start_pos)
                #print(f'var_windows: {type(var_windows)}')
                number += 1
                # 在这里进行逐个比对
                diff_array, meta_array = compare_kmers_with_meta(ref_windows, var_windows)                
                results.append({
                    'chromosome_group': group_name,
                    'sequence_id': seq_id,
                    'diff_array': diff_array,
                    'meta_array': meta_array
                })
            #print(f"all_kmer:{all_kmer}")
            #diff_array_test = compare_kmers_all(all_kmer, number) #写一个新函数

            #print(f'results:{results}')
            return {
                'file_name': file_name,
                'results': results,
                'error': None
            }
        
        else: # 对含poly_ins的组进行单独处理，是一定有seq2的
            #print(f'含poly_ins，group为{group_name}')
            
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
            
            # ref_windows：[{'pos': 906668, 'kmer': 'T'}, {'pos': 906669, 'kmer': '-'}]
            ref_windows = kmer_window_meta(reference_seq, start_pos)
            #print(f"ref_windows: {ref_windows}")
            #这边我加一个储存所有kmer的字典, 键：ref_windows, var_windows1,2,3,4...
            all_kmer = {} 
            all_kmer['ref_windows'] = ref_windows
            number = 1
            for seq_id, sequence in sequences:
                if seq_id == reference_id:
                    continue  
                    
                var_windows = kmer_window_meta(sequence, start_pos)
                #print(f'var_windows: {type(var_windows)}')
                all_kmer[f'var_windows_{number}'] = var_windows
                number += 1
                # 在这里进行逐个比对,得改成凑一组全部比对才对
                diff_array, meta_array = compare_kmers_with_meta(ref_windows, var_windows)
                
                results.append({
                    'chromosome_group': group_name,
                    'sequence_id': seq_id,
                    'diff_array': diff_array,
                    'meta_array': meta_array
                })
            #print(f"all_kmer:{all_kmer}")
            #diff_array_test = compare_kmers_all(all_kmer, number) #写一个新函数

            #print(f'results:{results}')
            return {
                'file_name': file_name,
                'results': results,
                'error': None
            }
    except Exception as e:
        logger.error(f"Error processing sequences in {file_name}: {str(e)}")
        return {'file_name': file_name, 'results': [], 'error': str(e)}

def compare_kmers_all(all_kmer: Dict, number: int):
    # number就是var的数量
    # all_kmer:{'ref_windows': [{'pos': 906670, 'kmer': 'A'}...], ...,
    # 'var_windows_2': [{'pos': 906670, 'kmer': 'A'}, {'pos': 906671, 'kmer': '-'}]
    diff_number = 0
    diff_array = []
    meta_array = []
    set_array = set() #用来检测是否有poly variant
    length = len(all_kmer['ref_windows'])
    for i in range(length):
        ref = all_kmer['ref_windows'][i]['kmer']
        set_array.add(ref) # 对集合
        for n in range(number):
            var = all_kmer[f'var_windows_{n + 1}'][i]['kmer']
            



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
            
            for future in tqdm(as_completed(futures), total=len(futures), desc="Generating kmers", bar_format="{desc}: {n_fmt}/{total_fmt} groups"):
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



from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
from tqdm import tqdm

def process_group(chrom, group):
    """
    处理每个chrom的组，返回处理结果。
    """
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
    retained_diff, retained_meta = retain_changed_columns_group_with_meta(diff_arrays, meta_arrays)
    if not retained_diff or not retained_meta or len(retained_diff) != len(sequence_ids):
        logger.warning(f"Invalid processing results for group {chrom}")
        return None
    
    # 处理结果，去掉 '-' 字符
    merged_results = []
    for i, seq_id in enumerate(sequence_ids):
        if i < len(retained_diff) and i < len(retained_meta):
            processed_meta = []
            for meta in retained_meta[i]:
                processed_meta.append({
                    'pos': meta['pos'],
                    'ref': meta['ref'],  # 可以移除 '-' 字符
                    'alt': meta['alt']   # 可以移除 '-' 字符
                })
            
            merged_results.append({
                'chromosome_group': chrom,
                'sequence_id': seq_id,
                'diff_array': retained_diff[i],
                'meta_array': processed_meta
            })
    
    return merged_results

def process_and_merge_results(input_csv: str, output_csv: str, threads: str):
    """
    Process and merge k-mer window results, removing '-' characters in the final results.
    """
    try:
        # 读取CSV文件
        df = pd.read_csv(input_csv)
        
        if df.empty:
            logger.warning("Input file is empty")
            return
        
        # 根据chromosome_group分组
        grouped = df.groupby('chromosome_group')
        
        # 计算进度条总任务数：所有组的数量
        total_groups = len(grouped)
        
        merged_results = []
        
        # 使用进度条
        with tqdm(total=total_groups, desc="Merging kmers", unit="group") as pbar:
            # 使用ThreadPoolExecutor来并行处理每个组
            with ThreadPoolExecutor(max_workers=threads) as executor:
                future_to_group = {executor.submit(process_group, chrom, group): chrom for chrom, group in grouped}
                
                # 遍历每个任务的结果
                for future in as_completed(future_to_group):
                    result = future.result()
                    if result:
                        merged_results.extend(result)
                    pbar.update(1)  # 每个任务完成后更新进度条
        
        if not merged_results:
            logger.warning("No valid merged results")
            return
        
        # 创建新的DataFrame并保存
        merged_df = pd.DataFrame(merged_results)
        merged_df.to_csv(output_csv, index=False)
        logger.info(f"Merged results saved to: {output_csv}")
        
    except Exception as e:
        logger.error(f"Error during processing: {str(e)}")
        raise
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


def save_kmer_results_to_csv(results: Dict, output_file: str) -> None:
    """Save kmer comparison results to a CSV file."""
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
    meta_rows: List[List[Dict]]
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

    return retained_diff, retained_meta

