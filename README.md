# Pangenome Heritability Tool

这是一个用于处理泛基因组结构变异并生成VCF(v4.2)格式文件的Python工具。该工具通过执行序列比对、k-mer窗口分析，并将结果转换为VCF格式，帮助分析泛基因组中的结构变异。

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

## 软件输入

- `VCF`文件及其索引文件：仅包含结构变异（SVs），`.vcf`文件或`.vcf.gz`文件以及其索引文件
- `FASTA`文件：包含SVs对应染色体的FASTA序列文件（即`.fa`文件或`.fasta`文件）

## 软件主要输出
- `alignment_error_logs`文件夹：MAFFT运行的错误log会存放在该文件夹
- `alignment_results`文件夹：存储各个重叠组别的align结果，详情请看文末对文件夹结构的详细介绍
- `matrix_results`文件夹：存储各个重叠组别的矩阵输出结果
- `merged_rSVs.csv`：rSV的中间缓存，如`run-all`功能中断，可从该文件中途运行，不必重新执行`run-all`
- `rSV.vcf`文件：重叠的SV被细化成为rSV（refined SV）的vcf文件
- `rSV_meta.csv`文件：rSV具体的对应的pos，ref，alt细节
- `nrSV_meta.csv`文件：nrSV（non-overlapped rSV）, 即rSV中没有成功对齐，在给定阈值`cutoff`（默认0.9）的重叠度之下的片段对应的pos, ref, alt细节
- `nSV.vcf`文件：没有重叠的SV与nrSV的vcf文件
- `*.log`文件：输出的log文件，包括了运行时间，总处理的SV数量，INV数量，nSV数量，重叠SV数量，重叠率，总分组数量和rSV的总数信息

## 注意事项

**重要**：
- VCF和参考FASTA文件必须使用数字染色体标识符（例如：1、2、3表示染色体，如23表示X染色体，24表示Y染色体），且不应有任何前缀（如chr）或后缀。确保您的文件遵循此格式，以避免处理错误。
- 确保VCF文件的SV格式是标准化的：删除（deletion）应使用下面的VCF文件示例中的`sv1`标准格式，插入（insertion）应使用`sv2`，倒位（inversion）应使用`sv3`。
- 请确保`FORMAT`仅有`GT`，如不是，请使用外部工具进行提取。
- 您可以使用外部工具如`bcftools norm`进行标准化。
- 请尽量分染色体进行rSV的识别运行，以防止计算性能瓶颈。


### VCF文件示例：
VCF文件需要压缩使用，即后缀为`.vcf.gz`；确保Ref和Alt字段中的indels遵循标准格式，并且必须有索引文件，否则，请使用bcftools等软件生成索引文件，可能的代码：`bcftools index your_vcf_files.vcf.gz`
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
## 命令概览
该工具提供了以下命令：
- `process-vcf`：处理VCF文件，分组重叠变异，仅生成无重叠的VCF文件。
- `run-all`：在一个命令中运行整个工作流。（主命令）
- `make-meta`: 在有merged_rSVs.csv的情况下若运行中断，可由此继续`run-all`的步骤
## 快速使用示例
```bash
panherit process-vcf \
    --vcf test.vcf.gz \
    --ref test.fasta \
    --out nSV_test \
    --threads 10
```
选项：
- `--vcf`：输入包含结构变异的VCF文件
- `--ref`：参考基因组FASTA文件
- `--out`：处理后的变异和FASTA文件的输出目录
- `--threads`：并行线程数（默认：10）

```bash
panherit run-all \
    --vcf test.vcf.gz \
    --ref test.fasta \
    --cutoff 0.9 \
    --out output_directory \
    --threads 8
```
选项：
- `--vcf`：输入包含结构变异的VCF文件
- `--ref`：参考基因组FASTA文件
- `--cutoff`：低于该重叠度的片段不会被当作rSV（默认0.9）
- `--out`：处理后的变异和FASTA文件的输出目录
- `--threads`：并行线程数（默认：10）

```bash
panherit make-meta \
    --vcf test.vcf.gz \
    --ref test.fasta \
    --cutoff 0.9 \
    --out output_directory \
    --threads 8
```
同上
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
## 代码思路
### 代码主要达成效果
<p align="center">
<img src="https://github.com/user-attachments/assets/69c43992-af28-4669-a35e-198771986241" width="800">
</p>

即，将重叠的SV通过对齐与窗口判断等手段将其划分为更加细致的SV（rSV）
### [Step 1] Processing VCF and generating FASTA...
<p align="center">
<img src="https://github.com/user-attachments/assets/40ad6400-d640-4518-8775-644a6d93972f" width="500">
</p>

整体思路是，将`INV`扔到`nSV.vcf`里面，再利用`ref, alt, pos`等信息判断SV是否重叠，未重叠的给到`nSV.vcf`，储存到变量`single_group`中，使用`filter_vcf`函数进行vcf文件生成。重叠的储存到变量`multi_group`中，利用pos信息进行pre-align，结果用`generate_fasta_sequences`储存到`variants_pre_aligned.fasta`文件中供下一步使用。

### [Step 2] Running alignments...
<p align="center">
<img src="https://github.com/user-attachments/assets/46d21d26-2a96-42d3-ab11-13f13d8cc13c" width="600">
</p>

这一步的算法下，只有`group has insertions with the same POS（poly-ins）`的情况下，才会对其变异组进行切片和MAFFT align。您可以通过比对`input`和`aligned`来检查MAFFT的比对情况。不含`poly-ins`的情况会被直接保存为`aligned.fasta`。可能的文件夹结果：

<p align="center">
<img src="https://github.com/user-attachments/assets/3fdbbaa0-3274-4e7e-8d20-967420563934" width="300">
</p>

其中`Group`后先是Group对应的`chrom`，后是group的`number`，然后是group第一个variant的`pos`。即：
- 有`input_origin.fasta`后缀的说明该组别有`poly-ins`即多个同pos的insertion的情况。
- `input_sliced_X.fasta`指第`X`个切片输入，对应的`aligned_sliced_X.fasta`指对应的MAFFT的align结果
- `aligned.fasta`为最终结果
### [Step 3] Generating rSVs...
#### 总步骤
<p align="center">
<img src="https://github.com/user-attachments/assets/544c06a0-1415-4518-92d3-73485fcbbe06" width="900">
</p>

#### 普通情况（无`poly-ins`）

<p align="center">
  <img src="https://github.com/user-attachments/assets/d2827349-1ca6-4c04-a4d4-3adf087ff458" width="600">
</p>
<p align="center"><b>Figure 1:</b> 普通情况下的rSV生成的流程图，其中del = deletion, ins = insertion</p>


- 首先对对齐后的每个纵向的窗口进行扫描
- 对变异（如`del`和`ins`）和参考序列`Ref`进行比对，相同则在`diff`中标记为`0`，不同则标记为`1`
- 对相邻的列进行合并，得到`diff_array`
- 划分最终的`rSV_meta`，即rSV对应的`pos, ref, alt`

#### 特殊情况（有`poly-ins`，并且有对齐后的变异重叠度未达到阈值）
<p align="center">
  <img src="https://github.com/user-attachments/assets/293738ae-9093-42bf-8237-977208a66e2a" width="600">
</p>
<p align="center"><b>Figure 1:</b> 特殊情况下的rSV生成的流程图，其中del = deletion, ins = insertion</p>

- 不同点在于，在扫描出现`poly-alt`，即初步merge后仍出现多alt的情况，如果重叠度没有达到阈值，则会被丢入`nrSV_meta.csv`中，最终被送入`nSV.vcf`
### [Step 4] Converting to VCF format and generating matrices...
<p align="center">
  <img src="https://github.com/user-attachments/assets/c1bf08e0-8566-4f56-ac8e-5d087b3e79b0" width="1200">
</p>
<p align="center"><b>Figure 1:</b> Step 4 总流程图</p>

#### Generating rSV_meta.csv with meta information (positions, reference, alternate alleles)...
- 即提取`merged_rSVs.csv`的`meta`信息，对`pos, ref, alt`进行标准化的过程
<p align="center">
  <img src="https://github.com/user-attachments/assets/45d7e0f1-80df-46a1-80a1-32260e13ee6f" width="1200">
</p>
<p align="center"><b>Figure 1:</b> rSV_meta.csv的示例</p>

#### Generating matrices...
##### D matrices（SV与rSV的对应关系矩阵）
<p align="center">
  <img src="https://github.com/user-attachments/assets/911a5ca9-9f2a-4ef0-9b7e-941940a95bbe" width="600">
</p>
<p align="center"><b>Figure 1:</b> D-matrix生成示例</p>

- 即从每一组的`merged_results.csv`中合并后的`diff_array`摘下来
- 可以看到每个`SV`都可以当作`rSV`的线性组合
#### X matrices（SV与样本基因型的对应关系矩阵）
<p align="center">
  <img src="https://github.com/user-attachments/assets/5d460532-493b-43f9-9194-68501cbf12a1" width="800">
</p>
<p align="center"><b>Figure 1:</b> X-matrix生成示例（Ind.=Individual）</p>

- 从每一组对应的`D-matrix`找到`<input>.vcf`文件中对应的SV的样本基因型数据，并进行对应编码
- 注意：`./.`编码为`-9`仅为图例，实际代码编码为`-999`
#### T matrices（rSV与样本基因型的对应关系矩阵）
<p align="center">
  <img src="https://github.com/user-attachments/assets/acd87cf4-7d83-43fa-855b-06b5101ed3c2" width="900">
</p>
<p align="center"><b>Figure 1:</b> T-matrix生成示例（Ind.=Individual）</p>

- 显然的，`T-matrix`可以从同组的`D-matrix`和`X-matrix`运算得到，而`genotype(GT)`矩阵，则为`VCF(v4.2)`格式的GT矩阵

### Generating rSV.vcf...
<p align="center">
  <img src="https://github.com/user-attachments/assets/5f8517b2-1d0e-4239-91cd-cdc3de10e57d" width="900">
</p>
<p align="center"><b>Figure 1:</b> rSV.vcf生成示例（Ind.=Individual）</p>

- 用前文的`rSV_meta`和`GT-matrix`就可以生成了





## 文件夹结构

### 1. 主文件夹
- **a. merged_rSV.csv**：储存rsv正确的ref和alt，同一个rsv会有多个alt
- **b. ⭐ rSV.vcf**：最终输出的rsv的VCF文件
- **c. ⭐ rsv_meta.csv**：填充入VCF文件中的ID，pos，ref，alt的初始文件
- **d. nSV.vcf**：最终输出的nSV的VCF文件
- **e. nrSV_meta.csv**: 可查阅nrSV(non-overlapped rSV)的信息（ref，pos，alt）
- **f. variants_pre_aligned.fasta**：按POS进行预比对过后的FASTA文件汇总
- **g. X.log**：log信息，示例：
```bash
2025-03-11 23:21:25 - INFO - Total runtime: 3:47:19
2025-03-11 23:21:25 - INFO - Total variants: 102,882
2025-03-11 23:21:25 - INFO - INV count: 330
2025-03-11 23:21:25 - INFO - nSV count: 73,738
2025-03-11 23:21:25 - INFO - Overlapping SVs: 29,144
2025-03-11 23:21:25 - INFO - Overlap percentage: 28.33%
2025-03-11 23:21:25 - INFO - Total variant groups: 8,119
2025-03-11 23:21:25 - INFO - Final rSV count: 49,658
```

### 2. 子文件夹：alignment_results
- **a. Group_input_origin.fasta**：Group_"chrom"_"number"_"pos"，从variants_pre_aligned.fasta截取并简化的FASTA文件，作为比对的输入
- **b. Group_aligned.fasta**：上述文件经过比对后的结果文件，用于下一步k-mer的生成
- **c. Group_input_spliced.fasta**：同一位置的插入序列的切片
- **d. Group_aligned_spliced.fasta**：切片完成后，MAFFT软件比对的结果，对应的无切片后缀的aligned.fasta文件就是合并后的最终比对结果

### 3. 子文件夹：alignment_error_logs
- MAFFT比对的错误信息会生成在此文件夹中

### 4. 子文件夹：matrix_results
- **a. Group_D_matrix.csv**：Group_"chrom"_"number"_"pos"，D矩阵，rSV-SV
- **b. Group_T_matrix.csv**：对应的T矩阵，rSV-samples
- **c. Group_X_matrix.csv**：对应的X矩阵，SV-samples
