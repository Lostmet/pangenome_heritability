import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import List, Dict
from tqdm import tqdm

from ..exceptions import AlignmentError
from ..utils.logging_utils import get_logger
from ..variant_processing.vcf_parser import VariantGroup

logger = get_logger(__name__)

class AlignmentResult:
    def __init__(self, group_id: str, sequences: Dict[str, str]):
        self.group_id = group_id
        self.sequences = sequences
        self.reference = sequences.get('reference', '')
        
    @property
    def variant_count(self) -> int:
        return len(self.sequences) - 1  # Excluding reference
        
def run_muscle(input_fasta: Path, output_fasta: Path) -> None:
    """Run MUSCLE alignment on a single group"""
    try:
        result = subprocess.run(
            ["muscle5", "-align", str(input_fasta), "-output", str(output_fasta)],
            capture_output=True,
            text=True,
            check=True
        )
    except subprocess.CalledProcessError as e:
        raise AlignmentError(
            f"MUSCLE alignment failed: {e.stderr}\n"
            f"Input file: {input_fasta}"
        )

def align_variant_group(group: VariantGroup, 
                       config,
                       temp_dir: Path) -> AlignmentResult:
    """Align a single variant group"""
    input_fasta = temp_dir / f"group_{group.chrom}_{group.start}.fa"
    output_fasta = temp_dir / f"group_{group.chrom}_{group.start}_aligned.fa"
    
    # Write sequences to FASTA
    with open(input_fasta, 'w') as f:
        f.write(f">reference\n{group.get_reference_sequence()}\n")
        for i, var in enumerate(group.variants):
            f.write(f">variant_{i}\n{var.get_sequence()}\n")
            
    # Run alignment
    run_muscle(input_fasta, output_fasta)
    
    # Parse results
    sequences = {}
    with open(output_fasta) as f:
        current_id = None
        current_seq = []
        for line in f:
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id:
            sequences[current_id] = ''.join(current_seq)
            
    return AlignmentResult(
        group_id=f"{group.chrom}_{group.start}",
        sequences=sequences
    )

def run_alignments(config, variant_groups: List[VariantGroup]) -> List[AlignmentResult]:
    """Run alignments for all variant groups"""
    temp_dir = Path(config.output_dir) / "temp_alignments"
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    results = []
    with ProcessPoolExecutor(max_workers=config.threads) as executor:
        futures = [
            executor.submit(align_variant_group, group, config, temp_dir)
            for group in variant_groups
        ]
        
        for future in tqdm(futures, desc="Running alignments"):
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                logger.error(f"Alignment failed: {str(e)}")
                raise
                
    return results