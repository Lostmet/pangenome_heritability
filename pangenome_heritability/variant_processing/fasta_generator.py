import pysam
from tqdm import tqdm
import concurrent.futures
from ..config import Config
from typing import Dict, List, Tuple
from .vcf_parser import VariantGroup, Variant
from ..utils.logging_utils import get_logger


logger = get_logger(__name__)


def get_max_insertions(variants: List[Variant], start: int) -> Tuple[Dict[int, int], bool]:
    max_insertions = {}
    pos_set = []
    has_insertion = False  # 记录是否有多Insertion

    for variant in variants:
        rel_start = variant.start - start  # 计算相对位置
        if len(variant.ref) == 1 and len(variant.alt[0]) > 1:  # 仅对插入变异计算
            max_insertions[rel_start] = max(max_insertions.get(rel_start, 0), len(variant.alt[0])) 
            pos_set.append(variant.start)
    if len(pos_set) != len(set(pos_set)):
        has_insertion = True

    return max_insertions, has_insertion

def adjust_reference_for_insertions(ref_seq: str, max_insertions: Dict[int, int]) -> str:
    """修改参考序列，插入最大 `Insertion` 需要的 `"-"`"""
    ref_list = list(ref_seq)
    ins = 0
    for pos, max_ins_length in sorted(max_insertions.items()):
        blank = "-" * (max_ins_length - 1)
        blank = list(blank)
        position = 0
        for items in blank:
            ref_list.insert(pos + ins + position, items)  # 插入"-"
            position += 1

        ins += max_ins_length - 1 ##看看能不能对齐，成了


    return "".join(ref_list)

def adjust_variants_for_insertions(ref_seq: str, max_insertions: Dict[int, int], variant: Variant, start: int):
    ref_list = list(ref_seq)
    ins = 0
    rel_start = variant.start - start
    rel_end = variant.end - start + 1
    poly_ins = {}

    if len(variant.ref) > 1 and len(variant.alt[0]) == 1:
        # deletion
        ref_list[rel_start:rel_end] = ['-'] * (len(variant.ref) - 1)

        for pos, max_ins_length in sorted(max_insertions.items()):
            blank = ['-'] * (max_ins_length - 1)
            ref_list[pos + ins:pos + ins] = blank
            ins += len(blank)

    elif len(variant.ref) == 1 and len(variant.alt[0]) > 1:
        # insertion
        for pos, max_ins_length in sorted(max_insertions.items()):
            if rel_start != pos:
                blank = ['-'] * (max_ins_length - 1)
                ref_list[pos + ins:pos + ins] = blank
                ins += len(blank)
            else:
                alt = list(variant.alt[0][1:])
                blank = ['-'] * (max_ins_length - len(variant.alt[0]))
                ref_list[pos + ins:pos + ins] = alt + blank
                if blank:
                    poly_ins["start"] = pos + ins
                    poly_ins["end"] = pos + ins + len(alt + blank)
                    poly_ins["pos"] = start + 1
                ins += len(alt + blank)

    else:
        # SNP / 替换
        ref_list[rel_start - 1:rel_end] = list(variant.alt[0])
        for pos, max_ins_length in sorted(max_insertions.items()):
            blank = ['-'] * (max_ins_length - 1)
            ref_list[pos + ins:pos + ins] = blank
            ins += len(blank)

    return ''.join(ref_list), poly_ins

def generate_fasta_sequences(config: Config, variant_groups: Dict[str, List[VariantGroup]], total_groups: int) -> Tuple[str, Dict[str, bool]]:
    
    output_fasta = f"{config.output_dir}/variants_pre_aligned.fasta"
    ref_genome = pysam.FastaFile(config.ref_fasta)
    has_insertion_dict = {}  # 记录Group是否包含Insertion
    poly_ins_list = []

    def process_group(chrom, i, group):

        start = group.start - 1
        end = group.end + 1
        ref_seq = ref_genome.fetch(chrom, start, end).upper()

        max_insertions, has_insertion = get_max_insertions(group.variants, start)
        has_insertion_dict[f"Group_{chrom}_{i}_{group.start}"] = has_insertion

        ref_seq_adjusted = adjust_reference_for_insertions(ref_seq, max_insertions)
        seq_length = len(ref_seq_adjusted)

        variant_results = []
        for variant in group.variants:
            var_seq, poly_ins = adjust_variants_for_insertions(ref_seq, max_insertions, variant, start)
            var_id = f"Variant_{chrom}_{i}_{variant.start}_{variant.end}"
            variant_results.append((var_id, var_seq, poly_ins))

        bp_to_process = (len(group.variants) + 1) * seq_length

        return chrom, i, group, ref_seq_adjusted, variant_results, bp_to_process

    try:
        with open(output_fasta, 'w') as fasta_out:
            with tqdm(total=total_groups, desc="Pre-aligning variant groups", unit="group") as pbar:  
                futures = []
                with concurrent.futures.ThreadPoolExecutor(max_workers=config.threads) as executor:
                    for chrom, groups in variant_groups.items():
                        for i, group in enumerate(groups, 1):
                            futures.append(executor.submit(process_group, chrom, i, group))

                    seen = set()  # 用来记录已经出现过的字典
                    group_bp_all = 0 # 用来记录所有组需要处理的bp和

                    for future in concurrent.futures.as_completed(futures):
                        chrom, i, group, ref_seq_adjusted, variant_results, group_bp = future.result()
                        group_bp_all += group_bp

                        # 写入参考序列
                        fasta_out.write(f">Group_{chrom}_{i}_{group.start}\n{ref_seq_adjusted}\n")

                        # 写入变异序列和poly_ins信息
                        for var_id, var_seq, poly_ins in variant_results:
                            fasta_out.write(f">{var_id}\n{var_seq}\n")

                            # 更新 poly_ins_list
                            if poly_ins:  # 过滤掉空字典
                                dict_frozenset = frozenset(poly_ins.items())
                                if dict_frozenset not in seen:
                                    poly_ins_list.append(poly_ins)  # 添加不重复的字典
                                    seen.add(dict_frozenset)  # 记录该字典已出现

                        pbar.update(1) 

        return output_fasta, has_insertion_dict, poly_ins_list # 返回Group是否包含Insertion信息，还有poly_ins信息

    finally:
        ref_genome.close()  # 关闭参考基因组文件
