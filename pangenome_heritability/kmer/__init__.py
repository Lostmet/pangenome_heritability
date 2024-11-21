"""K-mer window generation and comparison functionality."""

from .window_generator import (
    KmerWindow,
    generate_windows
)
from .comparison import (
    compare_windows,
    process_comparison_matrix
)

__all__ = [
    'KmerWindow',
    'generate_windows',
    'compare_windows',
    'process_comparison_matrix'
]