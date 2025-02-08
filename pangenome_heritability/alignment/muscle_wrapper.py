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


def run_mafft(input_fasta: Path, output_fasta: Path, log_dir: Path = None) -> None:
    """Run MAFFT alignment and convert output to uppercase."""
    try:
        command = ["mafft", "--thread", "20", str(input_fasta)]

        
        # 运行 MAFFT 并捕获输出
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )

        # 转换结果中的序列为大写
        with open(output_fasta, "w") as output_file:
            for line in result.stdout.splitlines():
                if line.startswith(">"):
                    output_file.write(line + "\n")  # 头信息不变
                else:
                    output_file.write(line.upper() + "\n")  # 转换序列部分为大写
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

def align_group(group_name: str, sequences: List[str], temp_dir: Path, log_dir: Path, has_insertion: bool) -> AlignmentResult:
    """Align a single group of sequences using MAFFT or just write input if no insertion is present."""
    input_fasta = temp_dir / f"{group_name}_input.fasta"
    output_fasta = temp_dir / f"{group_name}_aligned.fasta"

    # Write sequences to FASTA
    with open(input_fasta, "w") as f:
        for i, seq in enumerate(sequences):
            f.write(f">seq{i}\n{seq}\n")
    
    
    if has_insertion:

        # Run MAFFT alignment if there is insertion
        run_mafft(input_fasta, output_fasta, log_dir=log_dir)

    else:
        # Just copy the input sequence to the output, no alignment
        with open(output_fasta, "w") as f:
            for seq in sequences:
                f.write(f">seq{sequences.index(seq)}\n{seq}\n")

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



def run_alignments(config, fasta_file: str, has_insertion_dict: Dict[str, bool]) -> List[AlignmentResult]:

    """Run alignments for all groups in the given FASTA file using MAFFT when necessary."""
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
        futures = {}

        for group_name, sequences in group_sequences.items():
            has_insertion = has_insertion_dict.get(group_name, False)  # 直接使用传入的 `has_insertion`，而不是从 `dict` 里取值
            # 提交到线程池，进行 MAFFT 对齐或直接写入
            futures[executor.submit(align_group, group_name, sequences, temp_dir, log_dir, has_insertion)] = group_name

        
        # 处理异步任务的结果
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
