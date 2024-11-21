"""Genotype mapping and PLINK format conversion."""

from .genotype_mapper import GenotypeMapper
from .plink_converter import (
    PlinkFiles,
    convert_to_plink
)

__all__ = [
    'GenotypeMapper',
    'PlinkFiles',
    'convert_to_plink'
]