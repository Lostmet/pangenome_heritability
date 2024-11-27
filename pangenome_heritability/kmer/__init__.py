"""K-mer window generation and comparison functionality."""

from .window_generator import (
    kmer_window,
    process_sequences,
    compare_windows,
    process_fasta_files,
    save_kmer_results_to_csv
)
from .comparison import (
    process_comparison_results,
)

__all__ = [
    'kmer_window',
    'process_sequences',
    'compare_windows',
    'process_fasta_files',
    'save_kmer_results_to_csv',
    'process_comparison_results'
]