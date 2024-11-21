from typing import List, Dict
import pandas as pd
import numpy as np

def compare_windows(windows: Dict[str, List[KmerWindow]]) -> pd.DataFrame:
    """Compare k-mer windows to reference sequence"""
    results = []
    
    for group_id, group_windows in windows.items():
        # Get reference windows
        ref_windows = {w.position: w.kmer for w in group_windows 
                      if w.sequence_id.startswith('Group_')}
        
        # Compare variant windows to reference
        for window in group_windows:
            if not window.sequence_id.startswith('Group_'):
                comparison = int(window.kmer == ref_windows.get(window.position, ''))
                results.append({
                    'chromosome_group': group_id,
                    'sequence_id': window.sequence_id,
                    'position': window.position,
                    'comparison': comparison
                })
    
    return pd.DataFrame(results)

def process_comparison_matrix(comparisons: pd.DataFrame) -> pd.DataFrame:
    """Process comparison results into a matrix and remove redundant columns"""
    # Pivot to create comparison matrix
    matrix = comparisons.pivot(index=['chromosome_group', 'sequence_id'],
                             columns='position',
                             values='comparison')
    
    # Remove redundant columns
    unique_patterns = matrix.T.drop_duplicates()
    matrix = matrix[unique_patterns.index]
    
    # Reset index and format results
    matrix = matrix.reset_index()
    matrix['comparison'] = matrix.iloc[:, 2:].values.tolist()
    return matrix[['chromosome_group', 'sequence_id', 'comparison']]