from dataclasses import dataclass
from typing import List, Dict, Generator
import numpy as np
from tqdm import tqdm

from ..exceptions import WindowError
from ..utils.logging_utils import get_logger

logger = get_logger(__name__)

@dataclass
class KmerWindow:
    """Represents a k-mer window"""
    group_id: str
    position: int
    sequence: str
    is_reference: bool = False
    
    @property
    def size(self) -> int:
        return len(self.sequence)
        
    def matches(self, other: 'KmerWindow') -> bool:
        return self.sequence == other.sequence

class WindowComparison:
    """Represents a comparison between reference and variant windows"""
    def __init__(self, group_id: str, position: int, matches_ref: bool):
        self.group_id = group_id
        self.position = position
        self.matches_ref = matches_ref

def generate_windows(sequence: str, 
                    window_size: int) -> Generator[KmerWindow, None, None]:
    """Generate k-mer windows from a sequence"""
    if window_size < 1:
        raise WindowError("Window size must be positive")
        
    if not sequence:
        raise WindowError("Empty sequence")
        
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        if '-' * window_size != window:  # Skip all-gap windows
            yield KmerWindow(
                position=i,
                sequence=window
            )

def compare_sequences(reference: str,
                     variant: str,
                     window_size: int,
                     group_id: str) -> List[WindowComparison]:
    """Compare two sequences using sliding windows"""
    if len(reference) != len(variant):
        raise WindowError("Sequences must be of equal length")
        
    comparisons = []
    ref_windows = list(generate_windows(reference, window_size))
    var_windows = list(generate_windows(variant, window_size))
    
    for ref_win, var_win in zip(ref_windows, var_windows):
        comparisons.append(
            WindowComparison(
                group_id=group_id,
                position=ref_win.position,
                matches_ref=ref_win.matches(var_win)
            )
        )
        
    return comparisons

def process_alignments(alignments: List['AlignmentResult'],
                      config) -> List[WindowComparison]:
    """Process all alignments and generate window comparisons"""
    all_comparisons = []
    
    with tqdm(alignments, desc="Processing windows") as pbar:
        for alignment in pbar:
            try:
                reference = alignment.reference
                for var_id, variant in alignment.sequences.items():
                    if var_id != 'reference':
                        comparisons = compare_sequences(
                            reference=reference,
                            variant=variant,
                            window_size=config.window_size,
                            group_id=alignment.group_id
                        )
                        all_comparisons.extend(comparisons)
            except Exception as e:
                logger.error(
                    f"Error processing alignment {alignment.group_id}: {str(e)}"
                )
                raise WindowError(str(e))
                
    return all_comparisons