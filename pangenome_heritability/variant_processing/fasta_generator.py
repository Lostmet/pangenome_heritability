from typing import Dict, List, Tuple
from ..config import Config
from .vcf_parser import VariantGroup, Variant
import pysam

def generate_fasta_sequences(config: Config, variant_groups: Dict[str, List[VariantGroup]]) -> Tuple[str, Dict[str, bool]]:
    output_fasta = f"{config.output_dir}/variants_extended.fasta"
    ref_genome = pysam.FastaFile(config.ref_fasta)
    has_insertion_dict = {}  # 记录Group是否包含Insertion

    try:
        with open(output_fasta, 'w') as fasta_out:
            for chrom, groups in variant_groups.items():
                for i, group in enumerate(groups, 1):
                    start = group.start - 1
                    end = group.end

                    ref_seq = ref_genome.fetch(chrom, start, end).upper()  # 只提取 group 内的序列

                    # 计算Insertion信息
                    max_insertions, has_insertion = get_max_insertions(group.variants, start)
                    has_insertion_dict[f"Group_{chrom}_{i}"] = has_insertion  # 记录是否包含Insertion

                    # 调整参考序列
                    ref_seq = adjust_reference_for_insertions(ref_seq, max_insertions)

                    # 写入参考序列
                    fasta_out.write(f">Group_{chrom}_{i}\n{ref_seq}\n")

                    for variant in group.variants:
                        var_seq = generate_variant_sequence(ref_seq, variant, start, max_insertions)
                        var_id = f"Variant_{chrom}_{i}_{variant.start}_{variant.end}"
                        fasta_out.write(f">{var_id}\n{var_seq}\n")

        return output_fasta, has_insertion_dict  # 返回Group是否包含Insertion信息

    finally:
        ref_genome.close()


def get_max_insertions(variants: List[Variant], start: int) -> Tuple[Dict[int, int], bool]:
    max_insertions = {}
    has_insertion = False  # 记录是否有Insertion

    for variant in variants:
        rel_start = variant.start - start  # 计算相对位置
        if len(variant.ref) == 1 and len(variant.alt[0]) > 1:  # 仅对插入变异计算
            max_insertions[rel_start] = max(max_insertions.get(rel_start, 0), len(variant.alt[0]))
            has_insertion = True  # 发现Insertion

    return max_insertions, has_insertion

def adjust_reference_for_insertions(ref_seq: str, max_insertions: Dict[int, int]) -> str:
    """修改参考序列，插入最大 `Insertion` 需要的 `"-"`"""
    ref_list = list(ref_seq)

    for pos, max_ins_length in sorted(max_insertions.items(), reverse=True):
        ref_list.insert(pos + 1, "-" * (max_ins_length - 1))  # 插入 `REF-1` 个 `"-"`

    return "".join(ref_list)

def generate_variant_sequence(ref_seq: str, variant: Variant, start: int, max_insertions: Dict[int, int]) -> str:
    """Generate pre-aligned variant sequence"""
    rel_start = variant.start - start
    rel_end = variant.end - start + 1

    ref_part = ref_seq[:rel_start]
    post_part = ref_seq[rel_end:]

    if len(variant.ref) > 1 and len(variant.alt[0]) == 1:
        # **处理 Deletion（缺失）**
        alt_part = "-" * len(variant.ref)
    elif len(variant.ref) == 1 and len(variant.alt[0]) > 1:
        # **处理 Insertion（插入）**
        max_insert_size = max_insertions.get(rel_start, len(variant.alt[0]))
        alt_part = variant.alt[0] + "-" * (max_insert_size - len(variant.alt[0]))
    elif "<INV>" in variant.alt:
        # **处理倒位**
        alt_part = reverse_complement(ref_seq[rel_start:rel_end])
    else:
        # **SNP 直接替换**
        alt_part = variant.alt[0]

    aligned_seq = ref_part + alt_part + post_part


    return aligned_seq

def reverse_complement(seq: str) -> str:
    """Generate reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))
