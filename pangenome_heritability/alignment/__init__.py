"""Sequence alignment using MUSCLE."""

from .muscle_wrapper import (
    run_alignments,
)
from .alignment_processor import process_alignments

__all__ = [
    'run_alignments',
    'process_alignments'
]
