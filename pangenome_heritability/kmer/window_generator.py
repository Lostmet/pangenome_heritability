import os
import glob
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

def process_sequences(file_name: str, sequences: List[Tuple[str, str]], k: int = 4) -> Dict:
    """Process sequences and handle missing reference cases."""
    try:
        ref_sequence = None
        for seq_id, sequence in sequences:
            if seq_id == 'seq0':
                ref_sequence = sequence
                break
        
        if not ref_sequence:
            logger.warning(f"No reference sequence found in {file_name}")
            return {
                'file_name': file_name,
                'results': [{
                    'chromosome_group': file_name.replace('_aligned.fasta', '').replace('_input.fasta', ''),
                    'sequence_id': seq_id,
                    'comparison': [0]
                } for seq_id, _ in sequences if seq_id != 'seq0'],
                'error': None
            }
        
        ref_windows = kmer_window(ref_sequence, k)
        file_results = []
        
        for seq_id, sequence in sequences:
            if seq_id == 'seq0':
                continue
            var_windows = kmer_window(sequence, k)
            comparison = compare_windows(ref_windows, var_windows)
            file_results.append({
                'chromosome_group': file_name.replace('_aligned.fasta', '').replace('_input.fasta', ''),
                'sequence_id': seq_id,
                'comparison': comparison
            })
        
        return {'file_name': file_name, 'results': file_results, 'error': None}
    
    except Exception as e:
        logger.error(f"Error processing sequences in {file_name}: {str(e)}")
        return {'file_name': file_name, 'results': [], 'error': str(e)}

def process_fasta_files(
    directory: str,
    k: int = 4,
    max_workers: Optional[int] = None,
    output_file: str = None,
    error_log: str = None
) -> Dict:
    """Process FASTA files and generate comparison results."""
    try:
        if max_workers is None:
            max_workers = 10
        
        logger.info(f"Starting FASTA processing with {max_workers} workers")
        fasta_contents = read_fasta_files(directory)
        results = []
        errors = []
        
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(process_sequences, file_name, sequences, k): file_name
                for file_name, sequences in fasta_contents.items()
            }
            
            for future in tqdm(as_completed(futures), total=len(futures), desc="Processing files"):
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
    """Retain only columns where changes occur between consecutive positions."""
    if not rows:
        return []
    
    num_rows = len(rows)
    num_cols = len(rows[0])
    retained = [[row[0]] for row in rows]
    
    for col in range(1, num_cols):
        if any(rows[row][col] != rows[row][col - 1] for row in range(num_rows)):
            for row in range(num_rows):
                retained[row].append(rows[row][col])
    
    return retained


# Rest of the helper functions remain the same
def kmer_window(sequence: str, k: int = 4) -> List[str]:
    """Generate k-mer windows from a sequence."""
    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]

def compare_windows(ref_windows: List[str], var_windows: List[str]) -> List[int]:
    """Compare windows between reference and variant sequences."""
    return [1 if any(ref_window == var_window for ref_window in ref_windows) else 0 
            for var_window in var_windows]



def process_and_merge_results(input_file: str, output_file: str) -> None:
    """Process comparison results and merge similar patterns."""
    try:
        logger.info(f"Reading comparison results from: {input_file}")
        df = pd.read_csv(input_file)
        processed_data = []
        
        for group_name, group_data in df.groupby('chromosome_group'):
            comparisons = group_data['comparison'].apply(lambda x: eval(x) if isinstance(x, str) else x)
            matrix = np.array(comparisons.tolist())
            
            if len(group_data) == 1:
            
                non_zero_columns = np.any(matrix != 0, axis=0)
                matrix = matrix[:, non_zero_columns]
    
            
            if matrix.size == 0 or matrix.tolist() == [[0]]:
                processed_data.append(group_data.iloc[0])
            else:
                new_row = group_data.iloc[0].copy()
                new_row['comparison'] = str(matrix[0].tolist())
                processed_data.append(new_row)
            continue

            
            num_variants = matrix.shape[0]
            unique_columns = np.eye(num_variants, dtype=int)
            matrix = np.hstack((matrix, unique_columns))
            
            matrix = np.unique(matrix, axis=1)
            non_zero_columns = np.any(matrix != 0, axis=0)
            matrix = matrix[:, non_zero_columns]
            
            for i, row in enumerate(matrix):
                new_row = group_data.iloc[i].copy()
                new_row['comparison'] = str(row.tolist())
                processed_data.append(new_row)
        
        processed_df = pd.DataFrame(processed_data)
        processed_df.to_csv(output_file, index=False)
        logger.info(f"Processed results saved to: {output_file}")
    
    except Exception as e:
        logger.error(f"Error processing results: {str(e)}")
        raise

def process_chromosome_groups(input_csv: str, output_csv: str) -> None:
    """Process CSV results and retain changed columns for each chromosome group."""
    data = defaultdict(list)
    
    # Read CSV file
    with open(input_csv, newline='', encoding='utf-8') as csvfile:
        reader = pd.read_csv(csvfile)
        
        for _, row in reader.iterrows():
            chromosome_group = row['chromosome_group']
            data[chromosome_group].append({
                'sequence_id': row['sequence_id'],
                'comparison': eval(row['comparison'])  # Convert string to list
            })
    
    processed_data = []
    
    # Process each chromosome group
    for group, rows in data.items():
        comparisons = [row['comparison'] for row in rows]
        retained_comparisons = retain_changed_columns_group(comparisons)
        
        for i, row in enumerate(rows):
            processed_data.append({
                'chromosome_group': group,
                'sequence_id': row['sequence_id'],
                'comparison': retained_comparisons[i]
            })
    
    # Write processed data to output CSV
    processed_df = pd.DataFrame(processed_data)
    processed_df['comparison'] = processed_df['comparison'].apply(str)  # Save as string
    processed_df.to_csv(output_csv, index=False, encoding='utf-8')
    logger.info(f"Processed results saved to: {output_csv}")


def save_kmer_results_to_csv(results: Dict, output_file: str) -> None:
    """Save kmer comparison results to a CSV file.
    
    Args:
        results: Dictionary containing processed results and errors
        output_file: Path to output CSV file
    """
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