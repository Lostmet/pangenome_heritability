from typing import List, Dict
import pandas as pd
import numpy as np
import os

from ..utils.logging_utils import get_logger

logger = get_logger(__name__)

from .window_generator import rSV_window

def process_comparison_results(input_file: str, output_file: str):
    """Post-process comparison results for empty sequences."""
    try:
        df = pd.read_csv(input_file)

        
        mask = df['chromosome_group'].str.contains('_input.fasta', na=False)
        df.loc[mask, 'chromosome_group'] = df.loc[mask, 'chromosome_group'].str.replace('_input.fasta', '', regex=False)
        df.loc[mask, 'matches_ref'] = 0  

        
        df.to_csv(output_file, index=False)
        logger.info(f"Processed comparison results saved to {output_file}")
    except Exception as e:
        logger.error(f"Error processing comparison results: {str(e)}")
        raise
