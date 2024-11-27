"""K-mer window generation and comparison functionality."""

from .window_generator import (
    compare_windows,
    process_sequences,
    process_fasta_files,
    save_kmer_results_to_csv,
    kmer_window
)
from .comparison import (
    process_comparison_matrix
)

__all__ = [
    'kmer_window',
    'process_sequences',
    'compare_windows',
    'process_comparison_matrix',
    'process_fasta_files',
    'save_kmer_results_to_csv'
]