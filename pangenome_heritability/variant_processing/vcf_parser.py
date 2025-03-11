import pysam
from typing import Dict, List, Tuple
from dataclasses import dataclass
from tqdm import tqdm
import os
from ..exceptions import InputError
from ..utils.logging_utils import get_logger

logger = get_logger(__name__)

@dataclass
class Variant:
    chrom: str
    start: int
    end: int
    ref: str
    alt: List[str]
    samples: Dict[str, Tuple[int, int]]

class VariantGroup:
    def __init__(self, chrom: str):
        self.chrom = chrom
        self.variants = []
        self.start = float('inf')
        self.end = 0
        
    def add_variant(self, variant: Variant):
        self.variants.append(variant)
        self.start = min(self.start, variant.start)
        self.end = max(self.end, variant.end)
        
    def overlaps(self, variant: Variant) -> bool:
        return (variant.start <= self.end and 
                variant.end >= self.start)

def parse_chrom(chrom: str) -> int:
    """
    Convert the chromosome string to an integer that reflects a more natural
    ordering (e.g., 1 < 2 < ... < 10 < X < Y, etc.). 
    """
    # Remove common prefix like 'chr' if present
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    
    # Special handling for X, Y, M (mitochondrial)
    if chrom.upper() == "X":
        return 23
    elif chrom.upper() == "Y":
        return 24
    elif chrom.upper() in ["M", "MT"]:
        return 25
    
    # Try converting to integer, fallback if there's an error
    try:
        return int(chrom)
    except ValueError:
        return 999999

def process_variants(config) -> List[VariantGroup]:
    """
    Process a VCF file, sorting all variants by chromosome order 
    (using parse_chrom) and by their start coordinate, and then grouping 
    overlapping variants into VariantGroup objects.
    """
    try:
        vcf = pysam.VariantFile(config.vcf_file)
    except Exception as e:
        raise InputError(f"Failed to open VCF file: {str(e)}")
        
    all_variants = []
    single_group = [] 
    inv_group = [] # 加一个inv_group
    logger.info("Reading variants from VCF...")
    with tqdm(desc="Reading variants",unit='variants') as pbar:
        try:
            var_bp_all = 0
            var_bp_max = 0
            for record in vcf.fetch():
                #print(f'alts:{list(record.alts)}')
                if list(record.alts) == ['<INV>']:
                    #print("有inv嘛")
                    inv = Variant(
                        chrom=record.chrom,
                        start=record.pos,
                        end=record.pos + len(record.ref) - 1,
                        ref=record.ref,
                        alt=list(record.alts),
                        samples={s: record.samples[s]['GT'] for s in record.samples}
                    )
                    inv_group.append(inv) # 扔到新的group里面
                else:
                    variant = Variant(
                        chrom=record.chrom,
                        start=record.pos,
                        end=record.pos + len(record.ref) - 1,
                        ref=record.ref,
                        alt=list(record.alts),
                        samples={s: record.samples[s]['GT'] for s in record.samples}
                    )
                    var_bp = abs(len(variant.ref)-len(variant.alt[0]))
                    var_bp_max = max(var_bp_max, var_bp)
                    var_bp_all += var_bp # 暂时不管inv
                    all_variants.append(variant)
                pbar.update(1)
        except Exception as e:
            raise InputError(f"Error reading VCF file: {str(e)}")
        finally:
            vcf.close()
            
    logger.info("Sorting variants by chromosome (natural order) and start position...")
    # Sort by parsed chromosome value, then by start coordinate
    all_variants.sort(key=lambda v: (parse_chrom(v.chrom), v.start))
    # 筛选出所有vcf的variants后续准备进行group的判断
    logger.info("Grouping sorted variants...")
    variant_groups = []
    inv_groups = []
    current_group = None
    inv_current_group = None
    inv_count = len(inv_group)
    variant_count = len(all_variants) + inv_count
    # 处理inv_group
    for inv in inv_group:
        if inv_current_group is None:
            # First variant
            inv_current_group = VariantGroup(inv.chrom)
            inv_current_group.add_variant(inv)
        else:
            # If chromosome is different or they do not overlap, start a new group
            if (inv_current_group.chrom != inv.chrom or 
                not inv_current_group.overlaps(inv)):
                inv_groups.append(inv_current_group)
                inv_current_group = VariantGroup(inv.chrom)
                inv_current_group.add_variant(inv)
            else:
                # Still in the same group
                inv_current_group.add_variant(inv)
    click.echo(f"Excluding {inv_count} inversions") # 记录排除掉的inv的个数
    # Add the last group if it exists
    if inv_current_group is not None:
        inv_groups.append(inv_current_group)


    with tqdm(desc="Grouping variants", total=variant_count) as pbar:
        for variant in all_variants:
            if current_group is None:
                # First variant
                current_group = VariantGroup(variant.chrom)
                current_group.add_variant(variant)
            else:
                # If chromosome is different or they do not overlap, start a new group
                if (current_group.chrom != variant.chrom or 
                    not current_group.overlaps(variant)):
                    variant_groups.append(current_group)
                    current_group = VariantGroup(variant.chrom)
                    current_group.add_variant(variant)
                else:
                    # Still in the same group
                    current_group.add_variant(variant)

            pbar.update(1)
        
        # Add the last group if it exists
        if current_group is not None:
            variant_groups.append(current_group)
    
    logger.info(f"Processed {len(variant_groups)} variant groups")

    
    multi_group = []
    multi_var_bp = 0
    multi_var_bp_max = 0
    single_sv_count = 0
    multi_var_bp_all = 0
    variant_max = 0
    # 直接把inv_groups扔到single里面吃灰
    for i in range(len(inv_groups)):
        single_group.append(inv_groups[i])
        single_sv_count += 1
    ### multi_group1:{'chrom': '2', 'variants': [Variant(chrom='2', start=906670, end=906670, ref='A', alt=['ATATATATATATA'], samples={'SL001_SL001':
    for i in range(len(variant_groups)):
        if len(variant_groups[i].variants) == 1:
            single_group.append(variant_groups[i])
            single_sv_count += 1
        else:
            multi_group.append(variant_groups[i])
            for n in range(len(variant_groups[i].variants)):
                variant = variant_groups[i].variants[n]
                # print(f"这里：{variant.alt[0]}")
                ref_bp = len(variant.ref)
                alt_bp = len(variant.alt[0])
                multi_var_bp = abs(ref_bp - alt_bp)
                multi_var_bp_all += multi_var_bp
                if multi_var_bp_max < multi_var_bp:
                    multi_var_bp_max = multi_var_bp
                    variant_max = variant
                
    percentage_sv_overlapped = (1- single_sv_count/variant_count)*100 # 给到一个重叠的SV的百分比

    # print(f'single_group1:{single_group[0].__dict__}')
    # multi group，var未group的bp，single group，single_sv_count，multi后的var_bp，和max
    return multi_group, var_bp_all, var_bp_max, single_sv_count, multi_var_bp_all, multi_var_bp_max, percentage_sv_overlapped, variant_max, single_group, inv_count, variant_count

import pysam
import os
from tqdm import tqdm
import click
def filter_vcf(config, single_group):
    nvcf_name = "nSV.vcf"
    output_vcf = os.path.join(config.output_dir, nvcf_name)
    
    # 读取输入文件的header信息
    with pysam.VariantFile(config.vcf_file, "r") as vcf_in:
        header_text = str(vcf_in.header)
    
    # 手动写入header到输出文件
    with open(output_vcf, "w") as f_out:
        f_out.write(header_text)
        
        # 重新打开输入文件遍历记录
        with pysam.VariantFile(config.vcf_file, "r") as vcf_in:
            for i in tqdm(range(len(single_group)), desc="Generating nSV: ", unit='variants'):
                variant = single_group[i].variants[0]
                chrom, pos = variant.chrom, variant.start
                
                # 在染色体范围内搜索匹配位置的记录
                for rec in vcf_in.fetch(chrom, pos, pos + 1):
                    if rec.pos == pos:
                        # 转换为标准VCF行并写入
                        vcf_line = str(rec).strip()  # 去除前后空白
                        f_out.write(vcf_line + '\n')  # 确保换行符一致性
                        break
    click.echo(f"Non-overlapping SVs(nSVs) saved as {nvcf_name}")
    return nvcf_name


