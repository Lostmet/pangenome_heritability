"""K-mer window generation and comparison functionality."""

from .window_generator import (
process_fasta_files, process_chromosome_groups, process_and_merge_results, read_fasta_files
)
from .comparison import (
    process_comparison_results,
)

__all__ = [
    'kmer_window',
    'process_sequences',
    'compare_windows',
    'process_fasta_files',
    'process_comparison_results'
    'process_chromosome_groups',
    'process_and_merge_results',
    'read_fasta_files'
]