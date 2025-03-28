import re
import click
import subprocess
from tqdm import tqdm
from Bio import SeqIO
from pathlib import Path
from typing import List, Dict
from ..exceptions import AlignmentError
from ..utils.logging_utils import get_logger, log_tqdm_summary
from concurrent.futures import ProcessPoolExecutor, as_completed


logger = get_logger(__name__)


class AlignmentResult:
    def __init__(self, group_id: str, sequences: Dict[str, str]):
        self.group_id = group_id
        self.sequences = sequences
        self.reference = sequences.get('reference', '')

    @property
    def variant_count(self) -> int:
        return len(self.sequences) - 1  # Exclude reference


def read_fasta(fasta_file: str) -> Dict[str, List[str]]:
    """Read a FASTA file and group sequences by group header."""
    sequences = {}
    current_group = None
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id.startswith("Group"):
            current_group = record.id
            sequences[current_group] = [str(record.seq)]
        elif record.id.startswith("Variant") and current_group:
            sequences[current_group].append(str(record.seq))
    return sequences


def parse_fasta(fasta_file: str) -> Dict[str, str]:
    """Parse a FASTA file and return a dictionary of id → sequence."""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


def run_mafft(threads: int, input_fasta: Path, output_fasta: Path, log_dir: Path = None) -> None:
    """Run MAFFT and write alignment result to output file (converted to uppercase)."""
    try:
        command = ["mafft", "--thread", str(threads), str(input_fasta)]
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )

        with open(output_fasta, "w") as output_file:
            for line in result.stdout.splitlines():
                if line.startswith(">"):
                    output_file.write(line + "\n")
                else:
                    output_file.write(line.upper() + "\n")
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


def align_group(threads: str, group_name: str, sequences: List[str], temp_dir: Path, log_dir: Path, has_insertion: bool, poly_ins_list) -> AlignmentResult:
    """
    Align a single group using MAFFT.
    If no insertion, output as-is; otherwise, perform sliced alignments.
    """
    output_fasta = temp_dir / f"{group_name}_aligned.fasta"
    pos = int(re.search(r'(\d+)$', group_name).group(1))

    if has_insertion:
        input_fasta = temp_dir / f"{group_name}_input_origin.fasta"

        # Write original sequences to input FASTA
        with open(input_fasta, "w") as f:
            for i, seq in enumerate(sequences):
                f.write(f">seq{i}\n{seq}\n")
        origin_fasta = parse_fasta(input_fasta)

        # Filter relevant insertions for current group
        relevant_ins = [item for item in poly_ins_list if item['pos'] == pos]

        for i, ins in enumerate(relevant_ins):
            input_sliced = temp_dir / f"{group_name}_input_sliced_{i+1}.fasta"
            output_sliced = temp_dir / f"{group_name}_aligned_sliced_{i+1}.fasta"

            # Write sliced region to input
            with open(input_sliced, "w") as f:
                for j, seq in enumerate(sequences):
                    f.write(f">seq{j}\n{seq[ins['start']:ins['end']]}\n")

            # Run MAFFT on sliced region
            run_mafft(threads, input_sliced, output_sliced, log_dir=log_dir)
            sliced_fasta = parse_fasta(output_sliced)

            try:
                for key in origin_fasta:
                    if key in sliced_fasta:
                        ori_list = list(origin_fasta[key])
                        ori_list[ins['start']:ins['end']] = list(sliced_fasta[key])
                        origin_fasta[key] = "".join(ori_list)
            except:
                click.echo(f"Warning: Failed to insert sliced alignment in {group_name}")

        # Save updated aligned sequences
        with open(output_fasta, "w") as f:
            for i, seq_id in enumerate(origin_fasta):
                f.write(f">seq{i}\n{origin_fasta[seq_id]}\n")
    else:
        # No insertions: output directly without alignment
        with open(output_fasta, "w") as f:
            for i, seq in enumerate(sequences):
                f.write(f">seq{i}\n{seq}\n")

    # Parse aligned FASTA
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


def run_alignments(config, fasta_file: str, has_insertion_dict: Dict[str, bool], poly_ins_list) -> List[AlignmentResult]:
    """
    Run alignments for all groups using MAFFT, based on insertion status.
    """
    group_sequences = read_fasta(fasta_file)

    # Prepare directories
    temp_dir = Path(config.output_dir) / "alignment_results"
    log_dir = Path(config.output_dir) / "alignment_error_logs"
    temp_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    results = []

    with ProcessPoolExecutor(max_workers=config.threads) as executor, tqdm(
        total=len(group_sequences),
        desc="Processed Groups",
        unit='groups'
    ) as pbar:
        futures = {
            executor.submit(
                align_group,
                config.threads,
                group_name,
                sequences,
                temp_dir,
                log_dir,
                has_insertion_dict.get(group_name, False),
                poly_ins_list
            ): group_name
            for group_name, sequences in group_sequences.items()
        }

        for future in as_completed(futures):
            group_name = futures[future]
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                if log_dir:
                    with open(log_dir / f"{group_name}_error.log", "a") as log_file:
                        log_file.write(f"Alignment failed for {group_name}: {str(e)}\n")
            finally:
                pbar.update(1)
        pbar.close()
        log_tqdm_summary(pbar, logger)

    return results
