from Bio import SeqIO
from typing import Dict, List, Tuple

def process_alignments(alignment_dir: str) -> Dict[str, List[Tuple[str, str]]]:
    """Process MUSCLE alignment results"""
    alignments = {}
    
    for file_path in glob.glob(os.path.join(alignment_dir, "*_aligned.fasta")):
        sequences = []
        for record in SeqIO.parse(file_path, "fasta"):
            sequences.append((record.id, str(record.seq)))
        
        group_name = os.path.basename(file_path).replace("_aligned.fasta", "")
        alignments[group_name] = sequences
    
    return alignments