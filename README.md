# Pangenome Heritability Tool

A Python tool for processing pangenome structural variants and generating PLINK format files. This tool helps analyze structural variants in pangenomes by performing sequence alignment, k-mer window analysis, and converting results to VCF format.

## Quick Start

### Installation with Conda (Recommended)
```bash
# Create environment
conda create -n test_panherit python=3.8
conda activate test_panherit

# Install panherit （修改了一下git，git到我的fork上）
git clone https://github.com/Lostmet/pangenome_heritability.git
cd pangenome_heritability
pip install .

# Install MAFFT
conda install conda-forge::mafft

```


# Usage Guide

The Pangenome Heritability Tool provides several commands to process variants, perform alignments, and generate PLINK files. Each step can be run independently or as part of a pipeline.

## Command Overview

The tool provides four main commands: （现在只能用run-all，我还没有分开这些功能）
- `process-vcf`: Process VCF file and group overlapping variants
- `run-alignments`: Perform sequence alignments using MUSCLE
- `process-kmers`: Generate and analyze k-mer windows
- `convert-to-vcf`: Convert results to VCF format
- `run-all`: Run the entire workflow in one command

## Note

Important: The VCF and reference FASTA files must use numeric chromosome identifiers (e.g., 1, 2, 3 for chromosomes) without additional prefixes or suffixes. Ensure your files adhere to this convention to avoid processing errors.

Example of a VCF File Header: （VCF现在可以不解压了，需要保证indel的Ref和Alt是标准格式）

```##fileformat=VCFv4.2
##source=YourTool
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
1       12345   rs123   ACTA    A       50      PASS    .
2       67890   rs456   G       GAAC    99      PASS    .
```
Example of a FASTA File: （注意FASTA的输入的表头是>1这种类型的，不能是别的，如果是别的的话，需要自己调整）
```
>1
AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
>2
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG
```


## Quickly Usage
```bash
panherit run-all \
    --vcf test.vcf \
    --ref test.fasta \
    --out output_directory \
    --threads 4
```
Options:
- `--vcf`: Input VCF file containing structural variants
- `--ref`: Reference genome FASTA file
- `--out`: Output directory for processed variants and FASTA files
- `--window-size`: 现在不需要输入这个了，为固定值1
- `--threads`: Number of parallel threads (default: 1)


## Requirements

- Python 3.8+
- MAFFT V7.526
- External dependencies:
  - pandas
  - numpy
  - biopython
  - click
  - tqdm

# 文件夹结构

## 1. 主文件夹
- **a. Group_?_?_?_D_matrix.csv**: Group_"chrom"_"number"_"pos"，D矩阵，rSV-SV
- **b. Group_?_?_?_T_matrix.csv**: 对应的T矩阵，rSV-samples
- **c. Group_?_?_?_X_matrix.csv**: 对应的X矩阵, SV-samples
- **d. GT_matrix.csv**: T矩阵的合并（rsv.vcf的GT矩阵，格式为1/0类型）
- **e. ⭐ output_final_results.csv**: 储存rsv正确的ref和alt，同一个rsv会有多个alt
- **f. output.vcf**: 无GT矩阵的vcf
- **g. ⭐ pangenome_rsv.vcf**: 最终输出的rsv的vcf文件
- **h. (processed_) comparison_results.csv**: kmer比对的过程文件
- **i. rsv_meta.csv**: 填充入vcf文件中的ID，pos，ref，alt的初始文件（目前还不是真正对应的ref和alt）
- **j. T_matrix_abnormal_all.csv**: rSV-sample矩阵（T矩阵）中异常值占非零正常值的比例
- **k. T_matrix_abnormal.csv**: 具体group_name下的，每个T矩阵中非零值的数目，异常值的数目和比例
- **l. variants_extended.fasta**: 按照POS截取并分组的fasta文件汇总

## 2. 子文件夹：alignment_results
- **a. Group_?_?_?_input_origin.fasta**: Group_"chrom"_"number"_"pos"，从variants_extended.fasta截取并简化的fasta文件，作为align的输入，进行了预对齐，没有切片
- **b. Group_?_?_?_aligned.fasta**: 上一个文件经过比对后的结果文件，用于下一步kmer的生成
- **c. Group_?_?_?_input_spliced_?.fasta**: 顾名思义，就是同pos的insertion的切片的序号
- **d. Group_?_?_?_aligned_spliced_?.fasta**：切片完成后，mafft软件比对的结果，对应的无spliced后缀的aligned.fasta文件就是合并后的最终align结果

  
## 3. 子文件夹：logs
- 错误信息一部分会生成于此
