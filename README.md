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

# Install MUSCLE (Option 1: with sudo)
wget https://drive5.com/muscle5/muscle5.1.linux_intel64
chmod +x muscle5.1.linux_intel64
sudo mv muscle5.1.linux_intel64 /usr/local/bin/muscle5

# OR Install MUSCLE (Option 2: without sudo)
mkdir -p ~/local/bin
mv muscle5.1.linux_intel64 ~/local/bin/muscle5
echo 'export PATH="$HOME/local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

# Install package
git clone https://github.com/yourusername/pangenome-heritability.git
cd pangenome-heritability
conda install -c conda-forge pandas numpy biopython click tqdm
pip install -e .
```

### Basic Usage
```bash
panherit --vcf input.vcf --ref reference.fasta --out output_dir
```

## Installation Options

### Using Mamba (Faster Alternative)
```bash
# Install mamba
conda install -c conda-forge mamba

# Create environment
mamba create -n panherit python=3.8
mamba activate panherit

# Install dependencies
mamba install -c conda-forge pandas numpy biopython click tqdm
```

### Using environment.yml
```bash
# Create environment from file
mamba env create -f environment.yml
conda activate panherit
```

For detailed installation instructions, see [Installation Guide](docs/installation.md).

## Command Line Options

```bash
panherit --help

Options:
  --vcf TEXT         Input VCF file containing structural variants [required]
  --ref TEXT         Reference genome FASTA file [required]
  --out TEXT         Output directory [required]
  --threads INTEGER  Number of parallel processing threads [default: 1]
  --window-size INTEGER  K-mer window size [default: 4]
  --help            Show this message and exit
```

## Output Files

```
output_dir/
├── temp_alignments/        # MUSCLE alignment files
├── variants.bed           # PLINK binary file
├── variants.bim           # PLINK variant information
├── variants.fam           # PLINK sample information
└── pipeline.log          # Runtime log
```

## Python API Usage

```python
from pangenome_heritability.config import Config
from pangenome_heritability.variant_processing import process_variants
from pangenome_heritability.alignment import run_alignments
from pangenome_heritability.kmer import process_windows
from pangenome_heritability.genotype import convert_to_plink

# Configure and run pipeline
config = Config(
    vcf_file="input.vcf",
    ref_fasta="reference.fa",
    output_dir="results",
    threads=4,
    window_size=4
)

variants = process_variants(config)
alignments = run_alignments(config, variants)
windows = process_windows(config, alignments)
bfiles = convert_to_plink(config, windows)
```

## Documentation

- [Installation Guide](docs/installation.md)
- [Usage Guide](docs/usage.md)
- [API Reference](docs/api.md)

## Requirements

- Python 3.8+
- MUSCLE 5
- External dependencies:
  - pandas
  - numpy
  - biopython
  - click
  - tqdm

## Performance Tips

1. Use mamba for faster dependency installation
2. Adjust thread count based on available CPU cores
3. Monitor memory usage for large datasets
4. Ensure sufficient disk space for temporary files

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
- Email: [Your contact information]