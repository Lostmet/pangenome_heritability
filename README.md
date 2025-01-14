# Pangenome Heritability Tool

A Python tool for processing pangenome structural variants and generating PLINK format files. This tool helps analyze structural variants in pangenomes by performing sequence alignment, k-mer window analysis, and converting results to VCF format.

## Features

- Process VCF files containing structural variants
- Align variant sequences using MUSCLE
- K-mer based variant quantification
- Generate VCF format files (ped/map/bfile)
- Parallel processing support
- Progress tracking and logging

## Quick Start

### Installation with Conda (Recommended)
```bash
# Create environment
conda create -n panherit python=3.8
conda activate panherit

# Install panherit
git clone https://github.com/PeixiongYuan/pangenome_heritability.git
cd pangenome_heritability
pip install .

# Install MUSCLE
mkdir -p ~/local/bin
cd ~/local/bin
wget https://github.com/rcedgar/muscle/releases/download/v5.3/muscle-linux-x86.v5.3


chmod +x muscle-linux-x86.v5.3
mv muscle-linux-x86.v5.3 muscle

echo 'export PATH="$HOME/local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

# Install MAFFT
conda install conda-forge::mafft


```


# Usage Guide

The Pangenome Heritability Tool provides several commands to process variants, perform alignments, and generate PLINK files. Each step can be run independently or as part of a pipeline.

## Command Overview

The tool provides four main commands:
- `process-vcf`: Process VCF file and group overlapping variants
- `run-alignments`: Perform sequence alignments using MUSCLE
- `process-kmers`: Generate and analyze k-mer windows
- `convert-to-vcf`: Convert results to VCF format
- `run-all`: Run the entire workflow in one command

## Note

Important: The VCF and reference FASTA files must use numeric chromosome identifiers (e.g., 1, 2, 3 for chromosomes) without additional prefixes or suffixes. Ensure your files adhere to this convention to avoid processing errors.

Example of a VCF File Header:

```##fileformat=VCFv4.2
##source=YourTool
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
1       12345   rs123   A       T       50      PASS    .
2       67890   rs456   G       C       99      PASS    .
```
Example of a FASTA File:
```
>1
AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
>2
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG
```


## Quickly Usage
```bash
panherit run-all \
    --vcf input.vcf \
    --ref reference.fasta \
    --out output_directory \
    --window-size 4 \
    --threads 4
```
Options:
- `--vcf`: Input VCF file containing structural variants
- `--ref`: Reference genome FASTA file
- `--out`: Output directory for processed variants and FASTA files
- `--window-size`: Size of k-mer windows (default: 4)
- `--threads`: Number of parallel threads (default: 1)

## Detailed Usage
### Step 1: Process VCF File
```bash
# Process VCF and generate FASTA sequences
panherit process-vcf \
    --vcf input.vcf \
    --ref reference.fasta \
    --out output_directory
```
Options:
- `--vcf`: Input VCF file containing structural variants
- `--ref`: Reference genome FASTA file
- `--out`: Output directory for processed variants and FASTA files

### Step 2: Run Alignments
```bash
# Perform sequence alignments
panherit run-alignments \
    --grouped-variants output_directory/variants.fasta \
    --ref reference.fasta \
    --out alignments_directory \
    --threads 4
```
Options:
- `--grouped-variants`: FASTA file from previous step
- `--ref`: Reference genome FASTA file
- `--out`: Output directory for alignments
- `--threads`: Number of parallel threads (default: 1)

### Step 3: Process K-mer Windows
```bash
# Process k-mer windows
panherit process-kmers \
    --alignments temp_alignments \
    --grouped-variants output_directory/variants.fasta \
    --window-size 4 \
    --out kmers_directory
```
Options:
- `--alignments`: Directory containing alignment results
- `--grouped-variants`: FASTA file from previous step
- `--window-size`: Size of k-mer windows (default: 4)
- `--out`: Output directory for k-mer results

### Step 4: Convert to VCF
```bash
# Generate VCF files
panherit convert-to-vcf \
    --csv-file kmers_directory/output_final_results.csv \
    --vcf-file kmers_directory/test.vcf.gz \
    --grouped-variants output_directory/variants.fasta \ 
    --out vcf_files
```
Options:
- `--csv-file`: CSV file from previous step
- `--vcf-file`: VCF file from previous step
- `--grouped-variants`: FASTA file from previous step
- `--out`: Output directory for VCF files

## Output Files

Each step produces specific output files:

### Process VCF
```
output_directory/
├── variants.fasta      # Grouped variant sequences
└── variants.log       # Processing log file
```

### Alignments
```
/path/to/output/
├── temp_alignments/
│   ├── Group_2_59_input.fasta
│   ├── Group_2_59_aligned.fasta
├── error_logs/
│   ├── Group_2_59_input_error.log

```

### K-mer Windows
```
kmers_directory/
├── windows.csv       # K-mer window analysis
└── comparison.log    # Processing log file
```

### VCF Files
```
vcfiles/
├── pangenome.rsv.vcf
└──  pangenome.sv.vcf.gz
```

## Notes

- Each command will create its output directory if it doesn't exist
- Log files are generated for each step
- Use `--help` with any command for detailed options
- For large datasets, adjust thread count based on available resources

## Error Handling

If any step fails, the tool will:
1. Display an error message
2. Log the error details
3. Exit with a non-zero status code

Example error checking:
```bash
panherit process-vcf --vcf input.vcf --ref ref.fa --out output || {
    echo "VCF processing failed"
    exit 1
}
```


## Documentation

- [Installation Guide](docs/installation.md)
- [Usage Guide](docs/usage.md)
- [API Reference](docs/api.md)

## Requirements

- Python 3.8+
- MUSCLE 5
- MAFFT V7.526
- External dependencies:
  - pandas
  - numpy
  - biopython
  - click
  - tqdm

## Performance Tips


1. Adjust thread count based on available CPU cores
2. Ensure sufficient disk space for temporary files

## Troubleshooting

Common issues and solutions are documented in our [troubleshooting guide](docs/troubleshooting.md).

## Contributing

Contributions are welcome! Please read our [contributing guidelines](CONTRIBUTING.md).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this tool in your research, please cite:
```
[Citation information to be added]
```

## Contact

For questions and support:
- Open an issue on GitHub
- Email: yuanpeixiong@westlake.edu.cn
