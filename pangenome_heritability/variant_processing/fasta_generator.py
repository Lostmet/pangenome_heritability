from typing import Dict, List, Tuple
from ..config import Config
from .vcf_parser import VariantGroup, Variant
import pysam
import numpy as np

def generate_fasta_sequences(config: Config, variant_groups: Dict[str, List[VariantGroup]]) -> Tuple[str, Dict[str, bool]]:
    output_fasta = f"{config.output_dir}/variants_extended.fasta"
    ref_genome = pysam.FastaFile(config.ref_fasta)
    has_insertion_dict = {}  # 记录Group是否包含Insertion
    poly_ins_dict = []

    try:
        with open(output_fasta, 'w') as fasta_out:
            for chrom, groups in variant_groups.items():
                for i, group in enumerate(groups, 1):
                    start = group.start - 1
                    end = group.end

                    ref_seq = ref_genome.fetch(chrom, start, end).upper()  # 只提取 group 内的序列

                    # 计算Insertion信息
                    max_insertions, has_insertion = get_max_insertions(group.variants, start)
                    has_insertion_dict[f"Group_{chrom}_{i}_{group.start}"] = has_insertion  # 记录是否包含Insertion

                    # 调整参考序列 (更改为了adjusted)
                    ref_seq_adjusted = adjust_reference_for_insertions(ref_seq, max_insertions)

                    # 写入参考序列
                    fasta_out.write(f">Group_{chrom}_{i}_{group.start}\n{ref_seq_adjusted}\n")

                    for variant in group.variants:
                        var_seq, poly_ins = adjust_variants_for_insertions(ref_seq, max_insertions, variant, start)
                        var_id = f"Variant_{chrom}_{i}_{variant.start}_{variant.end}"
                        fasta_out.write(f">{var_id}\n{var_seq}\n")
                        poly_ins_dict.append(poly_ins)#创建一个[{start: end:}{start: end:}]，对每个都有
                    seen = set()  # 用来记录已经出现过的字典
                    poly_ins_list = []
                    for d in poly_ins_dict:
                        if d:  # 过滤掉空字典
                            # 将字典转换为 frozenset 后进行去重
                            dict_frozenset = frozenset(d.items())
                            if dict_frozenset not in seen:
                                poly_ins_list.append(d)  # 添加不重复的字典
                                seen.add(dict_frozenset)  # 记录该字典已出现

                    #print(poly_ins_list)                    

        return output_fasta, has_insertion_dict, poly_ins_list  # 返回Group是否包含Insertion信息，还有poly_ins信息

    finally:
        ref_genome.close()


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
        #print(pos_set)# 发现多insertion
              

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
        #print(f"位置：{pos + ins}")
        ins += max_ins_length - 1 ##看看能不能对齐，成了


    return "".join(ref_list)

def adjust_variants_for_insertions(ref_seq: str, max_insertions: Dict[int, int], variant, start) -> str:
    #原始参考序列
    ref_list = list(ref_seq)
    # 加一个处理好的ref
    ins = 0
    rel_start = variant.start - start
    rel_end = variant.end - start + 1

    poly_ins = {} #记录多polyINS
    #处理deletion
    if len(variant.ref) > 1 and len(variant.alt[0]) == 1:
        #ref_list[rel_start - 1] = variant.alt[0] # 防止别人不标准化
        ref_list[rel_start:rel_end - 1] = list((len(variant.ref)-1)*"-") #把对应的位置进行替换，我恨-1
        for pos, max_ins_length in sorted(max_insertions.items()):
            blank = "-" * (max_ins_length - 1)
            blank = list(blank)
            position = 0
            for items in blank:
                ref_list.insert(pos + ins + position, items)  # 插入"-"
                position += 1
            #print(f"位置：{pos + ins}")
            ins += max_ins_length - 1 ##看看能不能对齐，成了
    
    #处理insertion 
    elif len(variant.ref) == 1 and len(variant.alt[0]) > 1:
        #ref_list = list(ref_seq)
        #print(max_insertions)
        
        for pos, max_ins_length in sorted(max_insertions.items()):
            if rel_start != pos: #对max_insertion中，如果rel_start不是列表的pos，那就给占位符
                blank = "-" * (max_ins_length - 1)
                blank = list(blank)
                position = 0

                for items in blank:
                    ref_list.insert(pos + ins + position, items)  # 插入"-"
                    position += 1
                #print(f"位置：{pos + ins}")
                ins += max_ins_length - 1 ##看看能不能对齐，成了

            else: #真正修改的地方，加入了同pos多insertion的成功比对
                
                alt = list(variant.alt[0][1:])
                blank = "-" * (max_ins_length - len(variant.alt[0]))
                blank = list(blank)
                #print(blank)
                position = 0
                for items in alt:
                    ref_list.insert(pos + ins + position, items)
                    position += 1  # 插入 alt
                for items in blank:
                    ref_list.insert(pos + ins + position, items)
                    position += 1  # 插入 blank
                
                if blank != []: #得到同insertion的相对位置
                    poly_ins["start"] = pos + ins
                    poly_ins["end"] = pos + ins + position

                ins += max_ins_length - 1 ##看看能不能对齐

    # INV处理
    elif "<INV>" in variant.alt:
        ref_list[rel_start -1: rel_end] = reverse_complement(ref_list[rel_start -1: rel_end])
        for pos, max_ins_length in sorted(max_insertions.items()):
            blank = "-" * (max_ins_length - 1)
            blank = list(blank)
            position = 0
            for items in blank:
                ref_list.insert(pos + ins + position, items)  # 插入"-"
                position += 1
            #print(f"位置：{pos + ins}")
            ins += max_ins_length - 1 ##看看能不能对齐，成了
    # SNP和其他情况
    else:
        ref_list[rel_start:rel_end - 1] = list(variant.alt)
        for pos, max_ins_length in sorted(max_insertions.items()):
            blank = "-" * (max_ins_length - 1)
            blank = list(blank)
            position = 0
            for items in blank:
                ref_list.insert(pos + ins + position, items)  # 插入"-"
                position += 1
            #print(f"位置：{pos + ins}")
            ins += max_ins_length - 1 ##看看能不能对齐，成了       
        

    return "".join(ref_list), poly_ins      
         
#原始代码，不用了
def generate_variant_sequence(ref_seq: str, variant: Variant, start: int, max_insertions: Dict[int, int], ins_all: int) -> str:
    """Generate pre-aligned variant sequence"""
    rel_start = variant.start - start
    rel_end = variant.end - start + 1
    ins_len = 0

    if len(variant.ref) > 1 and len(variant.alt[0]) == 1:
        # **处理 Deletion（缺失）**
        ref_part = ref_seq[:(rel_start)]
        post_part = ref_seq[(rel_end + ins_all):]
       
        alt_part = "-" * len(variant.ref)

    elif len(variant.ref) == 1 and len(variant.alt[0]) > 1:
        # **处理 Insertion（插入）**
        ref_part = ref_seq[:(rel_start + ins_len)]
        post_part = ref_seq[(rel_end + ins_len):]        
        max_insert_size = max_insertions.get(rel_start, len(variant.alt[0]))
        alt_part = variant.alt[0] + "-" * (max_insert_size - len(variant.alt[0]))

        ref_part = ref_part[:- (len(alt_part) - 1)]

        ins_len += len(variant.alt[0]) - 1

    elif "<INV>" in variant.alt:
        # **处理倒位**
        alt_part = reverse_complement(ref_seq[(rel_start + ins_len):(rel_end + ins_len)])
    else:
        # **SNP 直接替换**
        alt_part = variant.alt[0]
        print(alt_part)
    aligned_seq = ref_part + alt_part + post_part

    return aligned_seq

def reverse_complement(seq: str) -> list:
    """Generate reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return list("".join(complement.get(base, base) for base in reversed(seq)))
