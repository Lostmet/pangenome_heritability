"""Sequence alignment using MUSCLE."""

from .muscle_wrapper import (
    run_alignments,
    align_group
)
from .alignment_processor import process_alignments

__all__ = [
    'run_alignments',
    'align_group',
    'process_alignments'
]
