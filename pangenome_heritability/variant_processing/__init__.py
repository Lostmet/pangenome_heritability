"""Variant processing and grouping functionality."""

from .vcf_parser import (
    VariantGroup,
    process_variants
)
from .fasta_generator import (
    generate_fasta_sequences,
    generate_variant_sequence,
    reverse_complement
)

__all__ = [
    'VariantGroup',
    'process_variants',
    'generate_fasta_sequences',
    'generate_variant_sequence',
    'reverse_complement'
]