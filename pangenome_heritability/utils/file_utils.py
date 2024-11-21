import os
import shutil
from pathlib import Path
from typing import Union, List

def ensure_directory(path: Union[str, Path]) -> Path:
    """
    Ensure a directory exists, creating it if necessary.
    
    Args:
        path: Directory path to ensure
        
    Returns:
        Path object of the ensured directory
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path

def cleanup_temp_files(directory: Union[str, Path], 
                      patterns: List[str] = None) -> None:
    """
    Clean up temporary files matching specified patterns.
    
    Args:
        directory: Directory containing temporary files
        patterns: List of glob patterns to match files for deletion
    """
    if patterns is None:
        patterns = ['*.tmp', '*.temp', '*_input.fasta']
        
    directory = Path(directory)
    for pattern in patterns:
        for file_path in directory.glob(pattern):
            try:
                if file_path.is_file():
                    file_path.unlink()
                elif file_path.is_dir():
                    shutil.rmtree(file_path)
            except Exception as e:
                logger.warning(f"Failed to remove {file_path}: {str(e)}")

def get_absolute_path(path: Union[str, Path]) -> Path:
    """
    Convert path to absolute path, expanding user and environment variables.
    
    Args:
        path: Path to convert
        
    Returns:
        Absolute Path object
    """
    return Path(os.path.expandvars(os.path.expanduser(str(path)))).resolve()

# pangenome_heritability/utils/logging_utils.py
import logging
import sys
from pathlib import Path
from typing import Union, Optional

def setup_logging(log_file: Optional[Union[str, Path]] = None,
                 level: int = logging.INFO) -> None:
    """
    Setup logging configuration.
    
    Args:
        log_file: Optional path to log file
        level: Logging level
    """
    # Create logger
    logger = logging.getLogger('pangenome_heritability')
    logger.setLevel(level)
    
    # Create formatters
    file_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    console_formatter = logging.Formatter(
        '%(message)s'
    )
    
    # Create console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    
    # Create file handler if log_file is specified
    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(str(log_file))
        file_handler.setLevel(level)
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

def get_logger(name: str = None) -> logging.Logger:
    """
    Get logger instance.