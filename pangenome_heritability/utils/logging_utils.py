# pangenome_heritability/utils/logging_utils.py
import logging
from pathlib import Path

def setup_logging(output_dir: Path) -> None:
    """Setup basic logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(output_dir / 'pipeline.log')
        ]
    )

def get_logger(name: str) -> logging.Logger:
    """Get a logger instance"""
    return logging.getLogger(name)

# pangenome_heritability/utils/file_utils.py
from pathlib import Path
import shutil
from typing import Union, List

def ensure_dir(path: Union[str, Path]) -> Path:
    """Ensure a directory exists"""
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path

def cleanup_temp_files(directory: Path, 
                      patterns: List[str] = None) -> None:
    """Clean up temporary files"""
    if patterns is None:
        patterns = ['*.tmp', '*.fa', '*.fasta']
    
    for pattern in patterns:
        for path in directory.glob(pattern):
            if path.is_file():
                path.unlink()
            elif path.is_dir():
                shutil.rmtree(path)