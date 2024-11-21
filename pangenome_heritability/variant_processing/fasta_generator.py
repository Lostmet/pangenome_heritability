from typing import Dict, List
from ..config import Config
from .vcf_parser import VariantGroup

def generate_fasta_sequences(config: Config, variant_groups: Dict[str, List[VariantGroup]]) -> str:
    """Generate FASTA file from variant groups"""
    output_fasta = f"{config.output_dir}/variants.fasta"
    ref_genome = pysam.FastaFile(config.ref_fasta)
    
    try:
        with open(output_fasta, 'w') as fasta_out:
            for chrom, groups in variant_groups.items():
                for i, group in enumerate(groups, 1):
                    # Write reference sequence
                    ref_seq = ref_genome.fetch(chrom, group.start-1, group.end)
                    fasta_out.write(f">Group_{chrom}_{i}\n{ref_seq}\n")
                    
                    # Write variant sequences
                    for variant in group.variants:
                        var_seq = generate_variant_sequence(ref_seq, variant, group.start)
                        var_id = f"Variant_{chrom}_{i}_{variant['start']}_{variant['end']}"
                        fasta_out.write(f">{var_id}\n{var_seq}\n")
        
        return output_fasta
    
    finally:
        ref_genome.close()

def generate_variant_sequence(ref_seq: str, variant: dict, group_start: int) -> str:
    """Generate variant sequence based on reference sequence"""
    rel_start = variant['start'] - group_start
    rel_end = variant['end'] - group_start + 1
    
    if "<INV>" in variant['alts']:
        # Handle inversion
        return (ref_seq[:rel_start] + 
                reverse_complement(ref_seq[rel_start:rel_end]) +
                ref_seq[rel_end:])
    else:
        # Handle other variant types
        return ref_seq[:rel_start] + variant['alts'][0] + ref_seq[rel_end:]

def reverse_complement(seq: str) -> str:
    """Generate reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))