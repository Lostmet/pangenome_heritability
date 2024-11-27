from typing import List, Dict
import pandas as pd
import numpy as np
import os

from .window_generator import KmerWindow

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