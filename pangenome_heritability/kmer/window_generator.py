import os
import glob
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict
from dataclasses import dataclass
from ..utils.logging_utils import get_logger
from ..alignment.muscle_wrapper import AlignmentResult

logger = get_logger(__name__)

@dataclass
class KmerWindow:
    """Represents a k-mer window"""
    group_id: str
    position: int
    sequence: str
    is_reference: bool = False

    def matches(self, other: 'KmerWindow') -> bool:
        return self.sequence == other.sequence


class WindowComparison:
    """Represents a comparison between reference and variant windows"""
    def __init__(self, group_id: str, position: int, matches_ref: bool):
        self.group_id = group_id
        self.position = position
        self.matches_ref = matches_ref


def read_fasta_files(directory: str) -> Dict[str, List[str]]:
    """Read aligned FASTA files and store sequences in a dictionary."""
    pattern = os.path.join(directory, 'Group_*_*_aligned.fasta')
    file_paths = glob.glob(pattern)
    file_paths.sort()
    logger.info(f"Found {len(file_paths)} files matching the pattern.")

    fasta_contents = {}
    for file_path in file_paths:
        sequences = []
        for record in SeqIO.parse(file_path, "fasta"):
            sequences.append((record.id, str(record.seq)))
        if sequences:
            fasta_contents[os.path.basename(file_path)] = sequences
        else:
            logger.warning(f"Empty file detected: {file_path}")

    return fasta_contents


def kmer_window(sequence: str, k: int) -> List[str]:
    """Generate k-mer windows from a sequence."""
    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]


def compare_windows(ref_windows: List[str], var_windows: List[str]) -> List[int]:
    """Compare windows between reference and variant."""
    comparison_results = []
    for var_window in var_windows:
        match = any(ref_window == var_window for ref_window in ref_windows)
        comparison_results.append(1 if match else 0)
    return comparison_results


def process_sequences(file_name: str, sequences: List[str], k: int) -> Dict:
    """Process a single file's sequences by generating k-mer windows."""
    ref_sequence = None
    for seq_id, sequence in sequences:
        if seq_id == 'seq0':  # Reference sequence
            ref_sequence = sequence
            break

    if not ref_sequence:
        logger.error(f"Reference sequence 'seq0' not found in {file_name}")
        return {
            'file_name': file_name,
            'results': [
                {
                    'chromosome_group': file_name.replace('_aligned.fasta', ''),
                    'sequence_id': seq_id,
                    'comparison': [0]
                } for seq_id, _ in sequences
            ],
            'error': "Reference sequence missing or invalid"
        }

    ref_windows = kmer_window(ref_sequence, k)
    file_results = []

    for seq_id, sequence in sequences:
        if seq_id == 'seq0':
            continue
        var_windows = kmer_window(sequence, k)
        comparison = compare_windows(ref_windows, var_windows)
        file_results.append({
            'chromosome_group': file_name.replace('_aligned.fasta', ''),
            'sequence_id': seq_id,
            'comparison': comparison
        })

    return {'file_name': file_name, 'results': file_results, 'error': None}


def process_fasta_files(directory: str, k: int, max_workers: int = None) -> List[Dict]:
    """Process all FASTA files in a directory using multi-threading."""
    fasta_contents = read_fasta_files(directory)
    results = []

    if max_workers is None:
        max_workers = max(1, os.cpu_count() - 1)

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_sequences, file_name, sequences, k): file_name
            for file_name, sequences in fasta_contents.items()
        }

        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing FASTA files"):
            file_result = future.result()
            if file_result['results']:
                results.extend(file_result['results'])
            elif file_result['error']:
                logger.error(f"Error processing {file_result['file_name']}: {file_result['error']}")

    return results


def save_kmer_results_to_csv(results: List[Dict], output_file: str):
    """Save K-mer comparison results to a CSV file."""
    data = []
    for result in results:
        for res in result['results']:
            for position, match in enumerate(res['comparison']):
                data.append({
                    'chromosome_group': res['chromosome_group'],
                    'sequence_id': res['sequence_id'],
                    'position': position,
                    'matches_ref': match
                })

    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)
    logger.info(f"Results saved to {output_file}")
