import pandas as pd
import numpy as np
from typing import Dict, Tuple

class GenotypeMapper:
    def __init__(self, vcf_file: str):
        self.vcf = pysam.VariantFile(vcf_file)
        self.samples = list(self.vcf.header.samples)
        
    def map_genotypes(self, comparisons: pd.DataFrame) -> Tuple[Dict, pd.DataFrame]:
        """Map window comparisons to genotypes"""
        variant_genotypes = self._load_vcf_genotypes()
        window_genotypes = {}
        
        for _, row in comparisons.iterrows():
            chrom = self._extract_chromosome(row['chromosome_group'])
            pos = self._extract_position(row['sequence_id'])
            comp_vector = eval(row['comparison'])
            
            # Initialize genotypes for this window
            window_gt = {sample: None for sample in self.samples}
            
            # Map comparisons to genotypes
            for i, comp in enumerate(comp_vector):
                if comp == 1:
                    key = (chrom, pos)
                    if key in variant_genotypes:
                        self._update_window_genotypes(
                            window_gt, 
                            variant_genotypes[key]
                        )
            
            window_genotypes[(row['chromosome_group'], row['sequence_id'])] = window_gt
            
        return window_genotypes
    
    def _load_vcf_genotypes(self) -> Dict:
        """Load genotypes from VCF file"""
        genotypes = {}
        for record in self.vcf.fetch():
            key = (record.chrom, record.pos)
            genotypes[key] = {
                sample: record.samples[sample]['GT']
                for sample in self.samples
            }
        return genotypes
    
    def _update_window_genotypes(self, 
                               window_gt: Dict, 
                               variant_gt: Dict) -> None:
        """Update window genotypes with variant genotypes"""
        for sample, gt in variant_gt.items():
            if window_gt[sample] is None:
                window_gt[sample] = gt
            elif gt is not None:
                window_gt[sample] = tuple(
                    max(a, b) for a, b in zip(window_gt[sample], gt)
                )
