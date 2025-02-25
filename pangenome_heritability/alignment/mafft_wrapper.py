import click
import subprocess
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed
from pathlib import Path
from typing import List, Dict
from tqdm import tqdm
from Bio import SeqIO
import re

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

def parse_fasta(fasta_file: str) -> Dict[str, List[str]]:
    # 新定义的，用来读取切片align的函数
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)

    #print(sequences)
    return sequences


def run_mafft(threads: int, input_fasta: Path, output_fasta: Path, log_dir: Path = None) -> None:
    """Run MAFFT alignment and convert output to uppercase."""
    try:
        command = ["mafft", "--thread", str(threads), str(input_fasta)]
        # print(f"MAFFT输入threads:{threads}")
        
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

def align_group(threads: str, group_name: str, sequences: List[str], temp_dir: Path, log_dir: Path, has_insertion: bool, poly_ins_list) -> AlignmentResult:
    """Align a single group of sequences using MAFFT or just write input if no insertion is present."""
    output_fasta = temp_dir / f"{group_name}_aligned.fasta"
    pos = int(re.search(r'(\d+)$', group_name).group(1))
    #print(f"align_group_pos = {pos}")
    #print(sequences)
    # Write sequences to FASTA    
    if has_insertion:
        input_fasta = temp_dir / f"{group_name}_input_origin.fasta"
        output_fasta = temp_dir / f"{group_name}_aligned.fasta"

        with open(input_fasta, "w") as f:
            for i, seq in enumerate(sequences):
                f.write(f">seq{i}\n{seq}\n") #写入初始input
        origin_fasta = parse_fasta(input_fasta)

        # 对于同pos多insertion的list进行剿灭
        #print(f"优化前ins_list:{poly_ins_list}, pos = {pos}, 数值类型：{type(pos)}")
        poly_ins_list = [item for item in poly_ins_list if item['pos'] == pos]
        #print(f"优化后的ins_list:{poly_ins_list}")
        for i, ins_start_end in enumerate(poly_ins_list):
            input_fasta_spliced = temp_dir / f"{group_name}_input_spliced_{i+1}.fasta"
            output_fasta_spliced = temp_dir / f"{group_name}_aligned_spliced_{i+1}.fasta"
            #print(ins_start_end)
            with open(input_fasta_spliced, "w") as f:
                for i, seq in enumerate(sequences):
                    f.write(f">seq{i}\n{seq[ins_start_end['start']:ins_start_end['end']]}\n") #进行切片
            run_mafft(threads, input_fasta_spliced, output_fasta_spliced, log_dir=log_dir)#切片部分进行align操作
            spliced_fasta = parse_fasta(output_fasta_spliced)
            #print(f"切片：{spliced_fasta}")
            try:
                for key in origin_fasta: # {seq0:xxxxxxx, seq1:xxxxxxx}
                    #print(f"seq_id:{key}")
                    if key in spliced_fasta:
                        spliced_list = list(spliced_fasta[key])
                        ori_list = list(origin_fasta[key])
                        ori_list[ins_start_end['start']:ins_start_end['end']] = spliced_list[:]

                        origin_fasta[key] = "".join(ori_list)
            except:
                click.echo(f"Warning: Cannot replace spliced align results in {group_name}") 
        #print(f"我真的输出了：{origin_fasta}")


        with open(output_fasta, "w") as f:#根据切片进行align更改操作
            for i, seq in enumerate(origin_fasta):
                f.write(f">seq{i}\n{origin_fasta[seq]}\n")#写入切片序列




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



def run_alignments(config, fasta_file: str, has_insertion_dict: Dict[str, bool], poly_ins_list) -> List[AlignmentResult]:

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
        unit='group'
    ) as pbar:
        futures = {}

        for group_name, sequences in group_sequences.items():
            has_insertion = has_insertion_dict.get(group_name, False)  # 直接使用传入的 `has_insertion`，而不是从 `dict` 里取值
            # 提交到线程池，进行 MAFFT 对齐或直接写入
            futures[executor.submit(align_group, config.threads, group_name, sequences, temp_dir, log_dir, has_insertion, poly_ins_list)] = group_name

        
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
