import os
import glob
import re
import ast
import numpy as np
import pandas as pd
from tqdm import tqdm
from typing import List
from itertools import combinations
from collections import defaultdict
from typing import List, Dict, Tuple
from ..utils.logging_utils import get_logger, log_tqdm_summary
from concurrent.futures import ProcessPoolExecutor, as_completed


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
        number = 1
        for seq_id, sequence in sequences:
            if seq_id == reference_id:
                continue      
            var_windows = rSV_window_meta(sequence, start_pos)
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
    genome_metadata: dict, 
    max_workers: int
) -> Dict:
    """
    Process FASTA files and generate comparison results, integrating genome metadata.
    """
    try:
        
        fasta_contents = read_fasta_files(directory)
        results = []
        errors = []
        
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(process_sequences, file_name, sequences, genome_metadata, has_insertion_dict): file_name
                for file_name, sequences in fasta_contents.items()
            }
            
            with tqdm(total=len(futures), desc="Generating rSVs", unit="groups") as pbar:
                for future in as_completed(futures):
                    result = future.result()
                    if result['error']:
                        errors.append(f"Error in {result['file_name']}: {result['error']}")
                    if result['results']:
                        results.extend(result['results'])
                    pbar.update(1)  
                pbar.close()
                log_tqdm_summary(pbar, logger)
        
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
    Process the groups for each chromosome and return the results.
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
    
    # Process the data using retain_changed_columns_group_with_meta
    retained_diff, retained_meta, nrSV_set = retain_changed_columns_group_with_meta(diff_arrays, meta_arrays, chrom, has_insertion_dict, cutoff)
    if not retained_diff or not retained_meta or len(retained_diff) != len(sequence_ids):
        logger.warning(f"Invalid processing results for group {chrom}")
        return None
    
    if nrSV_set == set():
        nrSV_merged_results = []
        nrSV_dict = {}
    else:
        nrSV_merged_results, nrSV_dict = nrSV_convert(chrom, nrSV_set)

    # Process results
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

    return merged_results, nrSV_merged_results, nrSV_dict



def nrSV_convert(chrom: str, nrSV_set: set) -> List[dict]:
    '''
    Convert a set of string-encoded SV records into structured metadata results.
    '''
    merged_results = []
    nrSV_step = 0  # A step counter for nrSV indexing
    nrSV_dict = {'group': chrom, 'SV': []}  # Dictionary for easy VCF sample lookup later

    # Convert each string to a dictionary
    sv_records = [ast.literal_eval(item) for item in nrSV_set]

    # Group records by their 'SV' field
    groups = defaultdict(list)
    for rec in sv_records:
        groups[rec['SV']].append(rec)
    
    # For each SV group, sort by position and build merged result
    for sv_value, records in sorted(groups.items()):
        nrSV_dict['SV'].append(sv_value)
        
        # Sort records by position
        records.sort(key=lambda x: x['pos'])
        pos = records[0]['pos']
        ref = records[0]['ref']

        merged = []  # List to store merged segments
        current = records[0].copy()

        if len(records) == 1:
            merged = records
        else:
            for rec in records[1:]:
                next_pos = current["pos"] + len(current["ref"])  # Expected next position for merging
                if rec["pos"] == next_pos:
                    # Merge insertion-type variants
                    current["ref"] += rec["ref"]
                    current["alt"] += rec["alt"]
                else:
                    # Not mergeable, save current and reset
                    merged.append(current)
                    current = rec.copy()
            # Append the last merged item
            merged.append(current)

        # Build meta_array by combining REF and ALT (removing '-' from insertions)
        meta_array = [{'pos': pos, 'ref': ref + rec['ref'].replace('-', ''), 'alt': ref + rec['alt'].replace('-', '')} for rec in merged[1:]]
        nrSV_count = len(meta_array)
        
        for i in range(nrSV_count):
            merged_results.append({
                'chromosome_group': chrom + f'_nrSV{i+1+nrSV_step}',
                'sequence_id': f"seq{sv_value}",  # Sequence ID derived from SV value
                'meta_array': meta_array[i]
            })
        
        nrSV_step += nrSV_count

    # Sort merged results by nrSV index in chromosome_group
    merged_results.sort(key=lambda x: int(x['chromosome_group'].split('_')[4][4:]))

    return merged_results, nrSV_dict

def process_and_merge_results(results, output_csv: str, threads: str,
                              has_insertion_dict: dict, cutoff: float, nrSV_csv: str):
    """
    Process and merge window results.
    """
    # has_insertion_dict: {'Group_2_1_1210539': True}
    try:
        df = pd.DataFrame(results['processed'])
    except Exception as e:
        logger.error(f"Error saving results to CSV: {str(e)}")
        raise
        
    try:
        if df.empty:
            logger.warning("Input file is empty")
            return
        
        # Group by chromosome_group
        grouped = df.groupby('chromosome_group')
        
        # Count total tasks for progress bar
        total_groups = len(grouped)
        
        merged_results = []
        nrSV_merged_results = []
        nrSV_list = []

        # Use progress bar
        with tqdm(total=total_groups, desc="Merging rSVs", unit="group") as pbar:
            # Use ProcessPoolExecutor to process each group in parallel
            with ProcessPoolExecutor(max_workers=threads) as executor:
                future_to_group = {
                    executor.submit(process_group, chrom, group, has_insertion_dict, cutoff): chrom
                    for chrom, group in grouped
                }

                # Iterate over completed tasks
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
                        logger.warning(f"Error in processing group {group}: {e}")
                        logger.warning(f"Tip: Please check whether the sequences in {group} are normalized properly.")
                    pbar.update(1)
            pbar.close()
            log_tqdm_summary(pbar, logger)

        if not merged_results:
            logger.warning("No valid merged results")
            return
        
        # Create merged DataFrame and save to CSV
        merged_df = pd.DataFrame(merged_results)
        merged_df.to_csv(output_csv, index=False)
        logger.info(f"Merged results of rSV saved to: {output_csv}")

        # Save nrSV merged results
        nrSV_merged_df = pd.DataFrame(nrSV_merged_results)
        nrSV_merged_df.to_csv(nrSV_csv, index=False)
        logger.info(f"Merged results of nrSV saved to: {nrSV_csv}")

    except Exception as e:
        logger.error(f"Error during processing: {str(e)}")
        raise

    return nrSV_list

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

def retain_changed_columns_group_with_meta(
    rows: List[List[int]], 
    meta_rows: List[List[Dict]],
    chrom: str,
    has_insertion_dict: Dict,
    cutoff: float
) -> Tuple[List[List[int]], List[List[Dict]]]:
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

    num_rows = len(rows)
    num_cols = len(rows[0])

    retained_diff = [[rows[i][0]] for i in range(num_rows)]
    retained_meta = [[meta_rows[i][0]] for i in range(num_rows)]
    current_col = 0

    for col in range(1, num_cols):
        is_same = True
        for row_i in range(num_rows):
            if rows[row_i][col] != rows[row_i][col-1]:
                is_same = False
                break
        
        if is_same:
            for row_i in range(num_rows):
                curr_meta = meta_rows[row_i][col]
                prev_meta = retained_meta[row_i][-1]
                merged_meta = {
                    'pos': prev_meta['pos'],
                    'ref': prev_meta['ref'] + curr_meta['ref'][-1],
                    'alt': prev_meta['alt'] + curr_meta['alt'][-1]
                }
                retained_meta[row_i][-1] = merged_meta
        else:
            for row_i in range(num_rows):
                retained_diff[row_i].append(rows[row_i][col])
                retained_meta[row_i].append(meta_rows[row_i][col])
            current_col += 1

    nrSV_set = set()
    new_rSV_count = 0

    if has_insertion_dict[chrom]:
        sv_count = len(retained_diff)
        rSV_count = len(retained_diff[0])
        rSV_index_small = 0
        poly_pos = 0

        for i in range(1, rSV_count):
            rSV_index = i + rSV_index_small
            alt_list = []
            alt_set = set()
           
            for j in range(sv_count):
                if retained_diff[j][rSV_index] == 1 and retained_meta[j][rSV_index]['alt'][0] != "-":
                    alt_origin = retained_meta[j][rSV_index]['alt']
                    alt_set.add(alt_origin)
                    alt_list.append(alt_origin)
            
            poly_alt_count = len(alt_set)
                
            if poly_alt_count <= 1:
                for j in range(sv_count):
                    retained_meta[j][rSV_index]['pos'] += poly_pos
            else:
                alt_set_filtered = filter_similar_sequences(alt_list, cutoff)
                if alt_set_filtered == set():
                    poly_pos -= len(retained_meta[0][rSV_index]['ref'])
                    for j in range(sv_count):
                        retained_diff[j].pop(rSV_index)
                        nrSV_meta = retained_meta[j][rSV_index]
                        nrSV_meta['SV'] = j+1
                        nrSV_set.add(str(nrSV_meta))
                        meta0 = retained_meta[j][0]
                        meta0['SV'] = j+1
                        nrSV_set.add(str(meta0))
                        retained_meta[j].pop(rSV_index)
                    rSV_index_small -= 1
                else:
                    new_rSV_count += len(alt_set_filtered) - 1
                    poly_alt = 0
                    alt_index_list = []
                    meta_list = []
                    alt_nrSV_index_set = set()

                    for j in range(sv_count):
                        alt = retained_meta[j][rSV_index]['alt']
                        meta = retained_meta[j][rSV_index]
                        if alt not in alt_set_filtered and alt[0] != "-":
                            meta['SV'] = j+1
                            nrSV_set.add(str(meta))
                            meta0 = retained_meta[j][0]
                            meta0['SV'] = j+1  
                            nrSV_set.add(str(meta0))                          
                            retained_diff[j][rSV_index] = 0
                            alt_nrSV_index_set.add(j)

                    for alt_filtered in alt_set_filtered:
                        found = False
                        for j in range(sv_count):
                            if retained_meta[j][rSV_index]['alt'] == alt_filtered and not found:
                                alt_index_list.append(j)
                                found = True
                                meta = retained_meta[j][rSV_index]
                                break
                        meta['pos'] += poly_pos
                        meta_list.append(meta)        

                    _count_dict = {}
                    for index, alt_index in enumerate(alt_index_list):
                        meta = meta_list[index]
                        alt_to_compare = retained_meta[alt_index][rSV_index + poly_alt]['alt']
                        for j in range(sv_count):
                            if rSV_index + poly_alt < len(retained_meta[j]):
                                alt = retained_meta[j][rSV_index + poly_alt]['alt']
                            else:
                                alt = None
                            alt_del = retained_meta[j][rSV_index + poly_alt - index]['alt']
                            if alt_del[0] == "-" or j in alt_nrSV_index_set:
                                if f'SV{j}' not in _count_dict:
                                    _count_dict[f'SV{j}'] = j
                                    retained_meta[j][rSV_index + poly_alt]['pos'] = meta['pos']
                                elif _count_dict[f'SV{j}'] == j:
                                    retained_diff[j][rSV_index + poly_alt] = 0
                                    retained_meta[j][rSV_index + poly_alt] = meta  
                            elif alt != alt_to_compare:
                                retained_diff[j][rSV_index + poly_alt] = 0
                                retained_meta[j][rSV_index + poly_alt] = meta
                            else:                               
                                retained_meta[j][rSV_index + poly_alt] = meta

    return retained_diff, retained_meta, nrSV_set

def hamming_similarity(seq1, seq2):
    """Calculate Hamming similarity between two equal-length DNA sequences."""
    assert len(seq1) == len(seq2), "Sequences must be of equal length"
    return sum(a == b for a, b in zip(seq1, seq2)) / len(seq1)

def filter_similar_sequences(dna_list, threshold):
    """Retain DNA sequences with similarity ≥ threshold to at least one other sequence."""
    n = len(dna_list)
    similarity_matrix = np.zeros((n, n))

    for i, j in combinations(range(n), 2):
        similarity_matrix[i, j] = similarity_matrix[j, i] = hamming_similarity(dna_list[i], dna_list[j])

    max_similarities = similarity_matrix.max(axis=1)
    filtered_sequences = {dna_list[i] for i in range(n) if max_similarities[i] >= threshold}

    return filtered_sequences
