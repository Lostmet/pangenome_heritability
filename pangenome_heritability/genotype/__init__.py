"""Genotype mapping and PLINK format conversion."""

from .genotype_mapper import convert_to_plink_with_variants
from .plink_converter import (
    PlinkFiles,
    convert_to_plink
)

__all__ = [
    'convert_to_plink_with_variants',
    'PlinkFiles',
    'convert_to_plink'
]