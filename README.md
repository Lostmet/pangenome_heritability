# Pangenome Heritability Tool

🔬 **Pangenome Heritability Tool** 是一个用于处理 **Structured Variants (SVs)** 并生成 **VCF (v4.2) 格式**文件的 Python 工具。  
它通过 **序列比对、窗口扫描和合并**，生成 **refined SVs (rSVs)** 并转换为 VCF 格式，帮助下游分析。

---

## 🚀 快速开始

### 📌 使用 Conda 安装（推荐）

```bash
# 创建环境
conda create -n test_panherit python=3.8
conda activate test_panherit

# 安装 panherit
git clone https://github.com/Lostmet/pangenome_heritability.git
cd pangenome_heritability
pip install .

# 安装 MAFFT
conda install conda-forge::mafft
```

---

## 📂 软件输入

🔹 **输入文件要求：**
- **VCF** 文件及索引文件 (`.vcf.gz`与`.vcf.gz.tbi`或`.vcf.gz.csi`)
- **FASTA** 文件 (`.fa` 或 `.fasta`)，需包含 SVs 对应的染色体序列

---

## 📤 软件主要输出

| **文件/文件夹**        | **说明** |
|----------------------|---------|
| `alignment_error_logs/` | 存放 MAFFT 运行的错误日志 |
| `alignment_results/` | 存储各个重叠组别的比对结果 |
| `matrix_results/` | 存储各个重叠组别的矩阵输出 |
| `merged_rSVs.csv` | rSV 的中间缓存，运行中断时可用于恢复 |
| **`rSV.vcf`** | 细化后的 rSV 变异 VCF 文件（最终输出） |
| `rSV_meta.csv` | 详细的 rSV 变异信息（pos, ref, alt） |
| `nrSV_meta.csv` | 未达到重叠度阈值的 rSV 变异信息 |
| `nSV.vcf` | 无重叠 SV 与 nrSV 的 VCF 文件 |
| `*.log` | 运行日志，包括处理时间、SV 数量等信息 |

---

## ⚠ 重要注意事项

💡 **数据格式要求**
- **VCF 和参考 FASTA 文件必须使用数字染色体标识符**
  - 例如：`1`、`2`、`3`（染色体），`23`（X 染色体），`24`（Y 染色体）
  - **不能** 使用 `chr1`、`chrX` 等格式
- **SV 格式需标准化**：
  - **Deletion（缺失）** → `sv1` 示例
  - **Insertion（插入）** → `sv2` 示例
  - **Inversion（倒位）** → `sv3` 示例
- **FORMAT 字段仅包含 `GT`**
  - 若不符合，请使用 `bcftools norm` 进行标准化
- **建议分染色体运行 rSV 识别**
  - 避免计算性能瓶颈

---

## 📜 VCF 文件示例

⚠ **VCF 文件必须压缩 (`.vcf.gz`)，并有索引文件 (`.tbi`或`.csi`)**
```bash
bcftools index your_vcf_file.vcf.gz
```

示例 VCF：
```bash
##fileformat=VCFv4.2
##source=YourTool
#CHROM  POS  ID    REF     ALT     QUAL    FILTER  INFO   FORMAT Sample1  Sample2  Sample3  Sample4
1       1    sv1   ACTA    A       50      PASS    .        GT     1/1      1/0      0/0      ./.
1       5    sv2   G       GAAC    99      PASS    .        GT     0/0      1/0      0/0      0/0
1       6    sv3   GCTAG   <INV>   98      PASS    .        GT     ./.      0/0      1/1      1/1
```

---

## 📖 命令概览

| **命令** | **功能** |
|----------|---------|
| `process-vcf` | 处理 VCF 文件，分组重叠变异，仅生成无重叠 VCF |
| **`run-all`** | 一键执行完整流程（主命令） |
| `make-meta` | 运行中断后，恢复 `run-all` |

### **运行示例**

####  **仅处理 VCF**
如果您只想简单查看自己的VCF文件中的重叠情况，可以选择这一步
```bash
panherit process-vcf --vcf test.vcf.gz --ref test.fasta --out nSV_test --threads 10
```

####  **运行完整流程**
主功能
```bash
panherit run-all --vcf test.vcf.gz --ref test.fasta --cutoff 0.9 --out output_directory --threads 8
```

####  **恢复运行**
如果主功能运行时间过长导致中断，且生成了`merged_rSVs.csv`，您可以通过以下步骤恢复运行。<br>
**请注意：只需要将`run-all`改为`make-meta`，不要改变输出文件夹名**
```bash
panherit make-meta --vcf test.vcf.gz --ref test.fasta --cutoff 0.9 --out output_directory --threads 8
```

---

## 🛠 环境要求

- **Python 3.8+**
- **MAFFT v7.526**
- **依赖库**：
  - `pandas`
  - `numpy`
  - `biopython`
  - `click`
  - `tqdm`

---

## 🎯 代码思路

### 🏗 代码主要达成效果
<p align="center">
<img src="https://github.com/user-attachments/assets/69c43992-af28-4669-a35e-198771986241" width="800">
</p>
<p align="center"><b>Figure 1:</b> 代码主要达成效果示意图</p>

即，通过对齐与窗口判断等手段，将重叠的SV划分为更加细致的SV（rSV）。

---

## [Step 1] Processing VCF and Generating FASTA
<p align="center">
<img src="https://github.com/user-attachments/assets/40ad6400-d640-4518-8775-644a6d93972f" width="500">
</p>
<p align="center"><b>Figure 2:</b> Processing VCF and Generating FASTA 流程图</p>

**主要流程**：
- 将`INV`（倒位）变异加入 `nSV.vcf` 。
- 判断`SV`是否重叠：
  - **未重叠** → 存入 `single_group` 并用 `filter_vcf` 生成 VCF 文件。
  - **重叠** → 存入 `multi_group`，用 `pos` 进行 `pre-align`，然后存为 `variants_pre_aligned.fasta` 供下一步使用。

---

## [Step 2] Running Alignments
<p align="center">
<img src="https://github.com/user-attachments/assets/46d21d26-2a96-42d3-ab11-13f13d8cc13c" width="600">
</p>
<p align="center"><b>Figure 3:</b> Running Alignments 流程图</p>

**对齐逻辑**：
1. 仅当 **group 中包含相同 `POS` 的插入突变 (`poly-ins`)** 时，才会进行切片 (`sliced`) 和 MAFFT 对齐。
2. 结果文件说明：
   - `input_origin.fasta`：存在 `poly-ins` 的变异组。
   - `input_sliced_X.fasta`：第 `X` 个切片的输入文件。
   - `aligned_sliced_X.fasta`：对应 `MAFFT` 的对齐结果。
   - `aligned.fasta`：最终对齐结果。

<p align="center">
<img src="https://github.com/user-attachments/assets/3fdbbaa0-3274-4e7e-8d20-967420563934" width="300">
</p>
<p align="center"><b>Figure 4:</b> alignment_results文件夹示意图</p>

---

## [Step 3] Generating rSVs

### 总流程
<p align="center">
<img src="https://github.com/user-attachments/assets/544c06a0-1415-4518-92d3-73485fcbbe06" width="900">
</p>
<p align="center"><b>Figure 5:</b> Generating rSVs 总流程示意图</p>

### 情况 1：普通情况（无 poly-ins）
<p align="center">
<img src="https://github.com/user-attachments/assets/8817adcd-63d3-41e9-91bb-2990bdb6e07d" width="600">
</p>
<p align="center"><b>Figure 6:</b> 无 poly-ins 情况下的 rSV 生成流程</p>

1. 扫描对齐后的窗口。
2. 比对变异（`del`、`ins`）与参考序列：
   - 相同 → `diff` 标记 `0`
   - 不同 → `diff` 标记 `1`
3. 计算 `diff_array` 并合并相邻列，最终生成 `rSV_meta`。

### 情况 2：特殊情况（poly-ins 存在 & 变异重叠度未达阈值）
<p align="center">
<img src="https://github.com/user-attachments/assets/faa5ff4a-2b47-4a7c-93bd-e7991dd318c2" width="600">
</p>
<p align="center"><b>Figure 7:</b> poly-ins 存在时的 rSV 生成流程</p>

- 若 `poly-alt` 出现且重叠度不足阈值，则存入 `nrSV_meta.csv`，最终加入 `nSV.vcf`。

---

## [Step 4] Converting to VCF Format and Generating Matrices

### 总流程
<p align="center">
<img src="https://github.com/user-attachments/assets/c1bf08e0-8566-4f56-ac8e-5d087b3e79b0" width="1200">
</p>
<p align="center"><b>Figure 8:</b> Step 4 总流程示意图</p>

### 1. 生成 `rSV_meta.csv`
- 从 `merged_rSVs.csv` 提取 `pos`、`ref`、`alt` 信息并标准化。

<p align="center">
<img src="https://github.com/user-attachments/assets/45d7e0f1-80df-46a1-80a1-32260e13ee6f" width="1200">
</p>
<p align="center"><b>Figure 9:</b> rSV_meta.csv 示例</p>

### 2. 生成矩阵

#### **D 矩阵（SV 与 rSV 的对应关系矩阵）**
<p align="center">
<img src="https://github.com/user-attachments/assets/911a5ca9-9f2a-4ef0-9b7e-941940a95bbe" width="600">
</p>
<p align="center"><b>Figure 10:</b> D-matrix 生成示例</p>

- 从 `merged_results.csv` 提取 `diff_array` 。
- 每个 SV 可视作 rSV 的线性组合。

#### **X 矩阵（SV 与样本基因型的对应关系矩阵）**
<p align="center">
<img src="https://github.com/user-attachments/assets/5d460532-493b-43f9-9194-68501cbf12a1" width="800">
</p>
<p align="center"><b>Figure 11:</b> X-matrix 生成示例</p>

- 从 `D-matrix` 对应的 VCF 文件提取样本基因型 (`GT`) 并编码。
- `./.` 编码为 `-9`（示例），实际代码中为 `-999`。

#### **T 矩阵（rSV 与样本基因型的对应关系矩阵）**
<p align="center">
<img src="https://github.com/user-attachments/assets/8aedc9d7-fd57-4356-85dd-7c7c4a328f6e" width="900">
</p>
<p align="center"><b>Figure 12:</b> T-matrix 生成示例</p>

- `T-matrix = D-matrix × X-matrix`
- `GT-matrix` 为 VCF (v4.2) 格式的 GT 矩阵。

---

## 生成 `rSV.vcf`
<p align="center">
<img src="https://github.com/user-attachments/assets/5f8517b2-1d0e-4239-91cd-cdc3de10e57d" width="900">
</p>
<p align="center"><b>Figure 13:</b> rSV.vcf 生成示例</p>

- `rSV_meta + GT-matrix` 直接生成 `rSV.vcf` 文件。

---

## 📂 文件结构概览

### **1️⃣ 主目录**
| **文件/文件夹**        | **说明** |
|----------------------|---------|
| `merged_rSV.csv` | rSV 变异的中间文件 |
| **`rSV.vcf`** | 细化后的 rSV 变异（最终 VCF） |
| `rSV_meta.csv` | rSV 详细信息（pos, ref, alt） |
| `nSV.vcf` | 无重叠 SV 变异的 VCF 文件 |
| `nrSV_meta.csv` | 未达到重叠度的 rSV 信息 |
| `variants_pre_aligned.fasta` | 预比对的 FASTA 文件 |
| `*.log` | 运行日志 |

### **2️⃣ alignment_results 文件夹**
| **文件** | **说明** |
|---------|---------|
| `Group_input_origin.fasta` | 变异组初始 FASTA |
| `Group_aligned.fasta` | 变异组比对后的 FASTA |
| `Group_input_spliced.fasta` | 插入变异的切片 |
| `Group_aligned_spliced.fasta` | 切片后比对的结果 |

### **3️⃣ matrix_results 文件夹**
| **文件** | **说明** |
|---------|---------|
| `Group_D_matrix.csv` | D 矩阵（rSV-SV 关系） |
| `Group_T_matrix.csv` | T 矩阵（rSV-样本 关系） |
| `Group_X_matrix.csv` | X 矩阵（SV-样本 关系） |

---

## 🎯 总结
**Pangenome Heritability Tool** 提供了一整套 **SV 处理、比对和转换 VCF** 的工具，适用于大规模基因组分析。  
请确保格式正确，并根据需要选择 `process-vcf`、`run-all` 或 `make-meta` 进行处理。

---
