"""Pangenome-based heritability estimation pipeline."""

from .config import Config
from . import variant_processing
from . import alignment
from . import kmer
from . import genotype
from . import heritability
from . import utils

__version__ = '0.1.0'