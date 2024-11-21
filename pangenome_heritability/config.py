import os
from dataclasses import dataclass
from pathlib import Path

@dataclass
class Config:
    """Configuration for the pipeline"""
    vcf_file: str
    ref_fasta: str
    output_dir: str
    threads: int = 1
    window_size: int = 4
    
    def __post_init__(self):
        """Validate inputs and create output directory"""
        if not os.path.exists(self.vcf_file):
            raise FileNotFoundError(f"VCF file not found: {self.vcf_file}")
        if not os.path.exists(self.ref_fasta):
            raise FileNotFoundError(f"Reference FASTA not found: {self.ref_fasta}")
        os.makedirs(self.output_dir, exist_ok=True)

# pangenome_heritability/kmer/comparison.py
def process_windows(config: Config, alignments: Dict) -> pd.DataFrame:
    """Process alignments into k-mer windows and map to genotypes"""
    results = []
    for group_id, sequences in alignments.items():
        ref_seq = sequences[0]  # Reference sequence
        windows = generate_windows(ref_seq, config.window_size)
        
        for seq in sequences[1:]:  # Variant sequences
            comparisons = compare_to_reference(seq, windows)
            results.extend(comparisons)
    
    return pd.DataFrame(results)