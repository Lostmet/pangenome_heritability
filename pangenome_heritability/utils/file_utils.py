import pysam
import logging
import click
import re

def check_vcf_vs_fasta(vcf_path: str, ref_path: str):
    logger = logging.getLogger()

    # STEP 1: 提取 VCF 中实际使用的染色体 ID（遍历记录）
    try:
        with pysam.VariantFile(vcf_path) as vcf:
            vcf_chroms = set(rec.chrom for rec in vcf.fetch())
    except Exception as e:
        logger.error(f"Failed to open or parse VCF file: {e}")
        raise click.Abort()

    if not vcf_chroms:
        logger.warning("No variants found in VCF. Skipping chromosome consistency check.")
        return

    # STEP 2: 去除 'chr' 前缀
    first_chrom = next(iter(vcf_chroms))
    if first_chrom.startswith("chr"):
        logger.info("Detected 'chr' prefix in VCF chromosomes. Stripping it from all chromosome names.")
        vcf_chroms = {chrom[3:] if chrom.startswith("chr") else chrom for chrom in vcf_chroms}

    # STEP 3: 特殊染色体编码
    numeric = sorted([int(c) for c in vcf_chroms if re.fullmatch(r"\d+", c)])
    non_numeric = sorted([c for c in vcf_chroms if not re.fullmatch(r"\d+", c)])
    max_chr = max(numeric) if numeric else 0
    encoding = {name: max_chr + i + 1 for i, name in enumerate(non_numeric)}

    if encoding:
        logger.info("Special chromosome name encoding (for ordering):")
        for name, code in encoding.items():
            logger.info(f"  - {name} => {code}")

    # STEP 4: 提取参考序列中的染色体
    try:
        ref = pysam.FastaFile(ref_path)
        ref_chroms = set(ref.references)
    except Exception as e:
        logger.error(f"Failed to open reference FASTA file: {e}")
        raise click.Abort()

    # STEP 5: 比对是否缺失
    missing = sorted(chrom for chrom in vcf_chroms if chrom not in ref_chroms)
    if missing:
        logger.error("The following chromosomes are present in the VCF but missing from the reference FASTA:")
        for chrom in missing:
            logger.error(f"  - {chrom} (expected '> {chrom}' in FASTA)")

        logger.error("Chromosomes available in the reference FASTA:")
        for ref_name in sorted(ref_chroms):
            logger.error(f"  > {ref_name}")

        logger.error("Aborting due to inconsistent chromosome naming between VCF and reference FASTA.")
        raise click.Abort()

