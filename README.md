# Pangenome Heritability Tool

🔬 **Pangenome Heritability Tool** 是一个用于结构变异（Structural Variants, SVs）处理的 Python 工具，支持生成符合 VCF v4.2 格式的变异数据。  
该工具通过序列比对、窗口扫描和合并流程，构建 refined SVs（rSVs），主要聚焦于 `DEL`（缺失）和 `INS`（插入）类型。

由于 `INV`（倒位）在真实数据中占比极低，当前版本未处理该类型，以提升效率与实用性。例如，Zhou 等人在番茄泛基因组研究中报告 `INV` 仅占 SV 的 **0.32%**（Zhou et al., *Nature*, 2022, 606: 527–534）。


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
如果主功能运行时由于运行时间过长等原因被中断，且生成了`merged_rSVs.csv`，您可以通过以下步骤恢复运行。<br>
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
*如果您无需深入了解算法细节，只需确保文件格式正确，直接运行 `run-all` 即可完成所有步骤，并生成细分后的 `rSV.vcf` 文件。若希望进一步了解算法实现，请参考：[Algorithmic Logic of rSV Software](https://github.com/Lostmet/Algorithmic_Logic_of_rSV_Software)。*
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

log示例：
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

### **4️⃣ alignment_error_log 文件夹**

存放MAFFT运行错误log，正常为空文件夹

---

## 🎯 总结
**Pangenome Heritability Tool** 提供了一整套 **SV 处理、比对和转换 VCF** 的工具，适用于大规模基因组分析。  
请确保格式正确，并根据需要选择 `process-vcf`、`run-all` 或 `make-meta` 进行处理。

---
## 📕 参考文献
1.	Zhou, Y., et al., Graph pangenome captures missing heritability and empowers tomato breeding. Nature, 2022. 606(7914): p. 527-534.

---
