import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed
from pathlib import Path
from typing import List, Dict
from tqdm import tqdm
from Bio import SeqIO

from ..exceptions import AlignmentError
from ..utils.logging_utils import get_logger

logger = get_logger(__name__)

class AlignmentResult:
    def __init__(self, group_id: str, sequences: Dict[str, str]):
        self.group_id = group_id
        self.sequences = sequences
        self.reference = sequences.get('reference', '')
        
    @property
    def variant_count(self) -> int:
        return len(self.sequences) - 1  # Excluding reference


def read_fasta(fasta_file: str) -> Dict[str, List[str]]:
    """Read FASTA file and group sequences by group name."""
    sequences = {}
    current_group = None
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id.startswith("Group"):
            current_group = record.id
            sequences[current_group] = [str(record.seq)]
        elif record.id.startswith("Variant") and current_group:
            sequences[current_group].append(str(record.seq))
    return sequences


def run_muscle(input_fasta: Path, output_fasta: Path, use_super5: bool = False, log_dir: Path = None, timeout: int = 600) -> None:
    """Run MUSCLE alignment on a single group, optionally using Super5 method, with timeout detection."""
    try:
        if use_super5:
            command = ["muscle", "-super5", str(input_fasta), "-output", str(output_fasta)] 
        else:
            command = ["muscle", "-align", str(input_fasta), "-output", str(output_fasta)]
        
        # Run MUSCLE with a timeout
        subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True,
            timeout=timeout  # Set timeout in seconds
        )
    except subprocess.TimeoutExpired:
        # Handle timeout error
        if log_dir:
            log_dir.mkdir(parents=True, exist_ok=True)
            error_log = log_dir / f"{input_fasta.stem}_timeout.log"
            with open(error_log, "w") as log_file:
                log_file.write(f"MUSCLE alignment timed out for {input_fasta} after {timeout} seconds.\n")
        
        raise AlignmentError(
            f"MUSCLE alignment timed out for {input_fasta}. Error log saved at: {error_log if log_dir else 'Not logged'}"
        )
    except subprocess.CalledProcessError as e:
        if log_dir:
            log_dir.mkdir(parents=True, exist_ok=True)
            error_log = log_dir / f"{input_fasta.stem}_error.log"
            with open(error_log, "w") as log_file:
                log_file.write(f"MUSCLE alignment failed for {input_fasta}:\n")
                log_file.write(e.stderr)
        
        raise AlignmentError(
            f"MUSCLE alignment failed for {input_fasta}. Error log saved at: {error_log if log_dir else 'Not logged'}"
        )


def run_mafft(input_fasta: Path, output_fasta: Path, log_dir: Path = None) -> None:
    """Run MAFFT alignment as a fallback method."""
    try:
        command = ["mafft", "--thread", "20", str(input_fasta)]
        with open(output_fasta, "w") as output_file:
            subprocess.run(
                command,
                stdout=output_file,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
    except subprocess.CalledProcessError as e:
        if log_dir:
            log_dir.mkdir(parents=True, exist_ok=True)
            error_log = log_dir / f"{input_fasta.stem}_mafft_error.log"
            with open(error_log, "w") as log_file:
                log_file.write(f"MAFFT alignment failed for {input_fasta}:\n")
                log_file.write(e.stderr)
        
        raise AlignmentError(
            f"MAFFT alignment failed for {input_fasta}. Error log saved at: {error_log if log_dir else 'Not logged'}"
        )


def align_group(group_name: str, sequences: List[str], temp_dir: Path, log_dir: Path, muscle_timeout: int = 600) -> AlignmentResult:
    """Align a single group of sequences, switching to Super5 if normal alignment fails, and to MAFFT if Super5 fails."""
    input_fasta = temp_dir / f"{group_name}_input.fasta"
    output_fasta = temp_dir / f"{group_name}_aligned.fasta"
    
    # Write sequences to FASTA
    with open(input_fasta, "w") as f:
        for i, seq in enumerate(sequences):
            f.write(f">seq{i}\n{seq}\n")
    
    try:
        # Try aligning with the normal MUSCLE method
        run_muscle(input_fasta, output_fasta, log_dir=log_dir, timeout=muscle_timeout)
    except AlignmentError:
        # Log failure without printing
        if log_dir:
            log_dir.mkdir(parents=True, exist_ok=True)
            with open(log_dir / f"{group_name}_error.log", "a") as log_file:
                log_file.write(f"Normal MUSCLE alignment failed or timed out for {group_name}. Retrying with Super5 mode.\n")
        
        try:
            # Retry with Super5 method
            run_muscle(input_fasta, output_fasta, use_super5=True, log_dir=log_dir, timeout=muscle_timeout)
        except AlignmentError:
            # Log second failure without printing
            if log_dir:
                with open(log_dir / f"{group_name}_error.log", "a") as log_file:
                    log_file.write(f"Super5 MUSCLE alignment also failed or timed out for {group_name}. Switching to MAFFT.\n")
            # Use MAFFT as fallback if Super5 also fails
            run_mafft(input_fasta, output_fasta, log_dir=log_dir)
    
    # Parse results
    aligned_sequences = {}
    with open(output_fasta) as f:
        current_id = None
        current_seq = []
        for line in f:
            if line.startswith(">"):
                if current_id:
                    aligned_sequences[current_id] = "".join(current_seq)
                current_id = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id:
            aligned_sequences[current_id] = "".join(current_seq)
    
    return AlignmentResult(
        group_id=group_name,
        sequences=aligned_sequences
    )


def run_alignments(config, fasta_file: str) -> List[AlignmentResult]:
    """Run alignments for all groups in the given FASTA file"""
    # Parse FASTA file
    group_sequences = read_fasta(fasta_file)
    
    # Temporary directories for alignments and logs
    temp_dir = Path(config.output_dir) / "alignment_results"
    log_dir = Path(config.output_dir) / "logs"
    temp_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    
    results = []
    
    # Customize tqdm format for simpler progress
    with ProcessPoolExecutor(max_workers=config.threads) as executor, tqdm(
        total=len(group_sequences),
        desc="Processed Groups",
        bar_format="{desc}: {n}/{total} groups"
    ) as pbar:
        futures = {
            executor.submit(align_group, group_name, sequences, temp_dir, log_dir): group_name
            for group_name, sequences in group_sequences.items()
        }
        
        for future in as_completed(futures):
            group_name = futures[future]
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                # Log failure without printing
                if log_dir:
                    with open(log_dir / f"{group_name}_error.log", "a") as log_file:
                        log_file.write(f"Alignment failed for {group_name}: {str(e)}\n")
            finally:
                # Update progress bar
                pbar.update(1)
                
    return results
