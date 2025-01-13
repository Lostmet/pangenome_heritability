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
    
    >Group_2_2
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
      "Group_2_2": {
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


def kmer_window_meta(sequence: str, k: int, start_offset: int) -> list:
    """
    Given a sequence and its start coordinate in the genome (start_offset),
    return [{pos: ..., kmer: ...}, ...], where pos is the absolute coordinate.
    """
    windows = []
    for i in range(len(sequence) - k + 1):
        windows.append({
            'pos': start_offset + i,  # Calculate absolute coordinate
            'kmer': sequence[i:i + k]
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

    for i in range(length):
        ref = ref_windows[i]
        var = var_windows[i]
        
        
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


def process_sequences(file_name: str, sequences: List[Tuple[str, str]], genome_metadata: dict, k: int = 4) -> Dict:
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
    try:
       
        group_name = file_name.replace('_aligned.fasta', '').replace('_input.fasta', '')
        
       
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
            
            
            for i in range(len(ref_seq) - k + 1):
                ref_kmer = ref_seq[i:i+k]
                var_kmer = var_seq[i:i+k]
                
                
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
        
       
        ref_windows = kmer_window_meta(reference_seq, k, start_pos)
        
       
        for seq_id, sequence in sequences:
            if seq_id == reference_id:
                continue  
                
            
            var_windows = kmer_window_meta(sequence, k, start_pos)
            
            
            diff_array, meta_array = compare_kmers_with_meta(ref_windows, var_windows)
            
            
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
    genome_metadata: dict,  # Add genome_metadata parameter
    k: int = 4,
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
                executor.submit(process_sequences, file_name, sequences, genome_metadata, k): file_name
                for file_name, sequences in fasta_contents.items()
            }
            
            for future in tqdm(as_completed(futures), total=len(futures), desc="Processing files", bar_format="{desc}: {n_fmt}/{total_fmt} groups"):
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


def process_and_merge_results(input_csv: str, output_csv: str):
    """
    Process and merge k-mer window results, removing '-' characters in the final results.
    """
    try:
        # Read CSV file
        df = pd.read_csv(input_csv)
        
        if df.empty:
            logger.warning("Input file is empty")
            return
        
        # Group by chromosome_group
        grouped = df.groupby('chromosome_group')
        
        merged_results = []
        
        for chrom, group in grouped:
            try:
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
                    continue
                
                # Process data using retain_changed_columns_group_with_meta
                retained_diff, retained_meta = retain_changed_columns_group_with_meta(diff_arrays, meta_arrays)
                
                if not retained_diff or not retained_meta or len(retained_diff) != len(sequence_ids):
                    logger.warning(f"Invalid processing results for group {chrom}")
                    continue
                
                # Process results, removing '-' characters
                for i, seq_id in enumerate(sequence_ids):
                    if i < len(retained_diff) and i < len(retained_meta):
                        # Process '-' characters in meta_array
                        processed_meta = []
                        for meta in retained_meta[i]:
                            processed_meta.append({
                                'pos': meta['pos'],
                                'ref': meta['ref'].replace('-', ''),  # Remove '-' in ref
                                'alt': meta['alt'].replace('-', '')   # Remove '-' in alt
                            })
                        
                        merged_results.append({
                            'chromosome_group': chrom,
                            'sequence_id': seq_id,
                            'diff_array': retained_diff[i],
                            'meta_array': processed_meta
                        })
            
            except Exception as e:
                logger.error(f"Error processing group {chrom}: {str(e)}")
                continue
        
        if not merged_results:
            logger.warning("No valid merged results")
            return
            
        # Create new DataFrame and save
        merged_df = pd.DataFrame(merged_results)
        merged_df.to_csv(output_csv, index=False)
        logger.info(f"Merged results saved to: {output_csv}")
        
    except Exception as e:
        logger.error(f"Error during processing: {str(e)}")
        raise


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
            processed_data.append({
                'chromosome_group': chromosome_group,
                'sequence_id': row['sequence_id'],
                'diff_array': diff_array,  # Use original diff_array
                'meta_array': meta_array   # Use original meta_array
            })
    
    # Write processed CSV
    processed_df = pd.DataFrame(processed_data)
    # Convert arrays to strings for storage
    processed_df['diff_array'] = processed_df['diff_array'].apply(str)
    processed_df['meta_array'] = processed_df['meta_array'].apply(str)
    processed_df.to_csv(output_csv, index=False, encoding='utf-8')
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
        raise ValueError("Mismatch between the lengths of difference arrays and metadata arrays")

    num_rows = len(rows)
    num_cols = len(rows[0])

    # Initialize result lists
    retained_diff = [[rows[i][0]] for i in range(num_rows)]
    retained_meta = [[meta_rows[i][0]] for i in range(num_rows)]
    
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

# Example usage:
"""
# Input example:
diff_arrays = [
    [0, 0, 1, 1],
    [1, 1, 0, 0]
]

meta_arrays = [
    [
        {'pos': 100, 'ref': 'ATG', 'alt': 'ATG'},
        {'pos': 101, 'ref': 'TGC', 'alt': 'TGC'},
        {'pos': 102, 'ref': 'GC-', 'alt': 'GTA'},  # Contains deletion
        {'pos': 103, 'ref': 'CAT', 'alt': 'TAT'}
    ],
    [
        {'pos': 100, 'ref': 'ATG', 'alt': 'CTG'},
        {'pos': 101, 'ref': 'TGC', 'alt': 'TTC'},
        {'pos': 102, 'ref': 'GCA', 'alt': 'GCA'},
        {'pos': 103, 'ref': 'CAT', 'alt': 'CAT'}
    ]
]

# After processing:
retained_diff = [
    [0, 1],
    [1, 0]
]

retained_meta = [
    [
        {'pos': 100, 'ref': 'ATGC', 'alt': 'ATGC'},
        {'pos': 103, 'ref': 'CAT', 'alt': 'TAT'}  # Skipped sequence containing deletion
    ],
    [
        {'pos': 100, 'ref': 'ATGC', 'alt': 'CTTC'},
        {'pos': 102, 'ref': 'GCAT', 'alt': 'GCAT'}
    ]
]
"""
