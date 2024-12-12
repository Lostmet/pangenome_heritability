# Pangenome Heritability Tool

A Python tool for processing pangenome structural variants and generating PLINK format files. This tool helps analyze structural variants in pangenomes by performing sequence alignment, k-mer window analysis, and converting results to PLINK format.

## Features

- Process VCF files containing structural variants
- Align variant sequences using MUSCLE
- K-mer based variant quantification
- Generate PLINK format files (bed/bim/fam)
- Parallel processing support
- Progress tracking and logging

## Quick Start

### Installation with Conda (Recommended)
```bash
# Create environment
conda create -n panherit python=3.8
conda activate panherit

# Install MUSCLE
mkdir ~/local/bin
wget https://github.com/rcedgar/muscle/releases/download/v5.3/muscle-linux-x86.v5.3
chmod +x muscle-linux-x86.v5.3
mv muscle-linux-x86.v5.3 ~/local/bin/muscle
echo 'export PATH="$HOME/local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

# Install MAFFT
conda install conda-forge::mafft

# Install panherit
git clone https://github.com/PeixiongYuan/pangenome_heritability.git
cd pangenome_heritability
conda install -c conda-forge pandas numpy biopython click tqdm
pip install -e .
```


# Usage Guide

The Pangenome Heritability Tool provides several commands to process variants, perform alignments, and generate PLINK files. Each step can be run independently or as part of a pipeline.

## Command Overview

The tool provides four main commands:
- `process-vcf`: Process VCF file and group overlapping variants
- `run-alignments`: Perform sequence alignments using MUSCLE
- `process-kmers`: Generate and analyze k-mer windows
- `convert-to-plink`: Convert results to PLINK format

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
    --window-size 4 \
    --out kmers_directory
```
Options:
- `--alignments`: Directory containing alignment results
- `--window-size`: Size of k-mer windows (default: 4)
- `--out`: Output directory for k-mer results

### Step 4: Convert to PLINK
```bash
# Generate PLINK files
panherit convert-to-plink \
    --kmer-results kmers_directory \
    --grouped-variants output_directory/variants.fasta \ 
    --out plink_files
```
Options:
- `--kmer-results`: Directory containing k-mer analysis results
- `--grouped-variants`: FASTA file from previous step
- `--out`: Output directory for PLINK files

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

### PLINK Files
```
plink_files/
├── variants.bed      # Binary genotype file
├── variants.bim      # Variant information file
└── variants.fam      # Sample information file
```

## Example Pipeline

Complete pipeline example:
```bash
# 1. Process VCF
panherit process-vcf \
    --vcf input.vcf \
    --ref reference.fasta \
    --out step1_output

# 2. Run alignments
panherit run-alignments \
    --grouped-variants step1_output/variants.fasta \
    --ref reference.fasta \
    --out step2_output \
    --threads 4

# 3. Process k-mers
panherit process-kmers \
    --alignments step2_output \
    --window-size 4 \
    --out step3_output

# 4. Generate PLINK files
panherit convert-to-plink \
    --kmer-results step3_output \
    --out final_output
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
- PLINK 1.90
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
