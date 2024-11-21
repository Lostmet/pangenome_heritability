import pysam
from typing import Dict, List, Tuple
from dataclasses import dataclass
from tqdm import tqdm

from ..exceptions import InputError
from ..utils.logging_utils import get_logger

logger = get_logger(__name__)

@dataclass
class Variant:
    chrom: str
    start: int
    end: int
    ref: str
    alt: List[str]
    samples: Dict[str, Tuple[int, int]]

class VariantGroup:
    def __init__(self, chrom: str):
        self.chrom = chrom
        self.variants = []
        self.start = float('inf')
        self.end = 0
        
    def add_variant(self, variant: Variant):
        self.variants.append(variant)
        self.start = min(self.start, variant.start)
        self.end = max(self.end, variant.end)
        
    def overlaps(self, variant: Variant) -> bool:
        return (variant.start <= self.end and 
                variant.end >= self.start)

def process_variants(config) -> List[VariantGroup]:
    """Process VCF file and group overlapping variants"""
    try:
        vcf = pysam.VariantFile(config.vcf_file)
    except Exception as e:
        raise InputError(f"Failed to open VCF file: {str(e)}")
        
    variant_groups = []
    current_group = None
    
    with tqdm(desc="Processing variants") as pbar:
        try:
            for record in vcf.fetch():
                variant = Variant(
                    chrom=record.chrom,
                    start=record.pos,
                    end=record.pos + len(record.ref) - 1,
                    ref=record.ref,
                    alt=list(record.alts),
                    samples={s: record.samples[s]['GT'] for s in record.samples}
                )
                
                if current_group is None:
                    current_group = VariantGroup(variant.chrom)
                    current_group.add_variant(variant)
                elif (current_group.chrom != variant.chrom or 
                      not current_group.overlaps(variant)):
                    variant_groups.append(current_group)
                    current_group = VariantGroup(variant.chrom)
                    current_group.add_variant(variant)
                else:
                    current_group.add_variant(variant)
                    
                pbar.update(1)
                
            if current_group is not None:
                variant_groups.append(current_group)
                
        except Exception as e:
            raise InputError(f"Error processing VCF file: {str(e)}")
        finally:
            vcf.close()
            
    logger.info(f"Processed {len(variant_groups)} variant groups")
    return variant_groups