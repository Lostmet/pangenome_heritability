# Pangenome Heritability Tool

这是一个用于处理泛基因组结构变异并生成PLINK格式文件的Python工具。该工具通过执行序列比对、k-mer窗口分析，并将结果转换为VCF格式，帮助分析泛基因组中的结构变异。

## 快速开始

### 使用Conda安装（推荐）
```bash
# 创建环境
conda create -n test_panherit python=3.8
conda activate test_panherit

# 安装panherit
git clone https://github.com/Lostmet/pangenome_heritability.git
cd pangenome_heritability
pip install .

# 安装MAFFT
conda install conda-forge::mafft
```

### 使用指南

泛基因组遗传工具提供了多个命令来处理变异、执行比对并生成PLINK文件。每个步骤可以单独执行，也可以作为工作流的一部分运行。
## 软件输入

- `VCF`文件及其索引文件：仅包含结构变异（SVs），`.vcf`文件或`.vcf.gz`文件以及其索引文件
- `FASTA`文件：包含SVs对应染色体的FASTA序列文件（即.fa文件或.fasta文件）

## 软件主要输出
### 1. 文件主要输出：
- `alignment`文件夹：存储各个重叠组别的align结果，详情请看文末对文件夹结构的详细介绍
- `matrix_results`文件夹：存储各个重叠组别的矩阵输出结果
- `rSV.vcf`文件：重叠的SV被细化成为rSV（refined SV）的vcf文件
- `rSV_meta.csv`文件：rSV具体的对应的pos，ref，alt细节
### 2. 命令行输出：
- `run-all`命令：运行全部
- 示例：（还没写）

## 命令概览

该工具提供了以下命令：
- `process-vcf`：处理VCF文件，分组重叠变异，生成无重叠的VCF文件。
- `run-all`：在一个命令中运行整个工作流。
- `make-meta`: 在有out_final_results.csv无rSV_meta.csv的情况下，继续run-all的步骤
- `make-vcf`：在有rSV_meta.csv的情况下，继续run-all的步骤

## 注意事项

**重要**：VCF和参考FASTA文件必须使用数字染色体标识符（例如：1、2、3表示染色体），且不应有任何前缀或后缀。确保您的文件遵循此格式，以避免处理错误。确保VCF文件的SV格式是标准化的：删除（deletion）应使用`sv1`，插入（insertion）应使用`sv2`，倒位（inversion）应使用`sv3`。您可以使用外部工具如`bcftools norm`进行标准化（本软件未安装此工具）。请尽量分染色体进行rSV的识别运行，以防止计算性能瓶颈。

### VCF文件头示例：
VCF文件现在可以直接使用，无需解压；确保Ref和Alt字段中的indels遵循标准格式，并且必须有索引文件，否则，请使用bcftools等软件生成索引文件，可能的代码：`bcftools index your_vcf_files.vcf.gz`
```bash
##fileformat=VCFv4.2
##source=YourTool
#CHROM  POS  ID    REF     ALT     QUAL    FILTER  INFO   FORMAT Sample1  Sample2  Sample3  Sample4
1       1    sv1   ACTA    A       50      PASS    .        GT     1/1      1/0      0/0      ./.
1       5    sv2   G       GAAC    99      PASS    .        GT     0/0      1/0      0/0      0/0
1       6    sv3   GCTAG   <INV>   98      PASS    .        GT     ./.      0/0      1/1      1/1
```

### FASTA文件示例：
注意：FASTA文件的头部必须是`>1`格式。如果是其他格式，需要自行调整。请自行将特殊染色体，如X、Y、MT编码为具体数字
```
>1
ACTAGGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
>2
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG
```

## 快速使用示例
```bash
panherit process-vcf \
    --vcf test.vcf.gz \
    --ref test.fasta \
    --out nSV_test \
    --threads 8
```

选项：
- `--vcf`：输入包含结构变异的VCF文件
- `--ref`：参考基因组FASTA文件
- `--out`：处理后的变异和FASTA文件的输出目录
- `--threads`：并行线程数（默认：8）

```bash
panherit run-all \
    --vcf test.vcf.gz \
    --ref test.fasta \
    --out output_directory \
    --threads 8
```

选项：
- `--vcf`：输入包含结构变异的VCF文件
- `--ref`：参考基因组FASTA文件
- `--out`：处理后的变异和FASTA文件的输出目录
- `--threads`：并行线程数（默认：8）

## 环境要求
- Python 3.8+
- MAFFT V7.526
- 外部依赖：
  - pandas
  - numpy
  - biopython
  - click
  - tqdm

---

## 文件夹结构

### 1. 主文件夹
- **a. Group_D_matrix.csv**：Group_"chrom"_"number"_"pos"，D矩阵，rSV-SV
- **b. Group_T_matrix.csv**：对应的T矩阵，rSV-samples
- **c. Group_X_matrix.csv**：对应的X矩阵，SV-samples
- **d. GT_matrix.csv**：合并后的T矩阵（rsv.vcf的GT矩阵，格式为1/0类型）
- **e. output_final_results.csv**：储存rsv正确的ref和alt，同一个rsv会有多个alt
- **f. output.vcf**：无GT矩阵的VCF文件
- **g. ⭐ pangenome_rSV.vcf**：最终输出的rsv的VCF文件
- **h. comparison_results.csv**：k-mer比对的过程文件
- **i. ⭐ rsv_meta.csv**：填充入VCF文件中的ID，pos，ref，alt的初始文件
- **j. T_matrix_abnormal_all.csv**：rSV-sample矩阵（T矩阵）中异常值占非零正常值的比例
- **k. T_matrix_abnormal.csv**：具体group_name下的，每个T矩阵中非零值的数目，异常值的数目和比例
- **l. variants_extended.fasta**：按POS截取并分组的FASTA文件汇总

### 2. 子文件夹：alignment_results
- **a. Group_input_origin.fasta**：Group_"chrom"_"number"_"pos"，从variants_extended.fasta截取并简化的FASTA文件，作为比对的输入
- **b. Group_aligned.fasta**：上述文件经过比对后的结果文件，用于下一步k-mer的生成
- **c. Group_input_spliced.fasta**：同一位置的插入序列的切片
- **d. Group_aligned_spliced.fasta**：切片完成后，MAFFT软件比对的结果，对应的无切片后缀的aligned.fasta文件就是合并后的最终比对结果

### 3. 子文件夹：logs
- 错误信息部分会生成在此文件夹中
