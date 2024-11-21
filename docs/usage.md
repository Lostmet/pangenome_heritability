
## Usage Guide

### Basic Usage
```bash
panherit --vcf input.vcf --ref reference.fasta --out output_dir
```

### Parameters
- `--vcf`: Input VCF file containing structural variants
- `--ref`: Reference genome FASTA file
- `--out`: Output directory
- `--threads`: Number of parallel processing threads (default: 1)
- `--window-size`: k-mer window size (default: 4)

### Output Structure
```
output_dir/
├── temp_alignments/        # MUSCLE alignment temporary files
├── variants.bed           # PLINK binary file
├── variants.bim           # PLINK variant information file
├── variants.fam           # PLINK sample information file
└── pipeline.log          # Runtime log
```

### Example
Run with sample data:
```bash
panherit \
    --vcf examples/sample.vcf.gz \
    --ref examples/reference.fa \
    --out results \
    --threads 4
```

### Important Notes
1. Input File Requirements:
   - VCF file must contain structural variant information
   - VCF file should be bgzipped and indexed
   - Reference genome must match VCF coordinate system

2. Memory Usage:
   - Memory consumption scales with input data size and window size
   - For large datasets, consider increasing available memory

3. Temporary Files:
   - Temporary files are generated in the output directory during processing
   - Automatic cleanup after successful completion

## Troubleshooting

### 1. MUSCLE Not Found
Error message:
```
Command 'muscle5' not found
```
Solution:
- Verify MUSCLE installation
- Ensure muscle5 is in system PATH

### 2. Memory Error
Error message:
```
MemoryError: Unable to allocate array
```
Solution:
- Reduce number of parallel threads
- Increase system available memory
- Consider batch processing for large datasets

### 3. VCF File Format Error
Error message:
```
InputError: Failed to open VCF file
```
Solution:
- Verify VCF file format
- Check file compression and indexing
- Validate VCF file integrity with bcftools

## Code Examples

### Python API Usage
```python
from pangenome_heritability.config import Config
from pangenome_heritability.variant_processing import process_variants
from pangenome_heritability.alignment import run_alignments
from pangenome_heritability.kmer import process_windows
from pangenome_heritability.genotype import convert_to_plink

# Configure parameters
config = Config(
    vcf_file="input.vcf",
    ref_fasta="reference.fa",
    output_dir="results",
    threads=4,
    window_size=4
)

# Run pipeline
variants = process_variants(config)
alignments = run_alignments(config, variants)
windows = process_windows(config, alignments)
bfiles = convert_to_plink(config, windows)
```

## Performance Optimization Tips

1. Parallel Processing:
   - Increase thread count (--threads) to improve processing speed
   - Recommended not to exceed CPU core count

2. Memory Management:
   - Window size affects memory usage, adjust as needed
   - Consider batch processing for large datasets

3. Disk Space:
   - Ensure sufficient space in output directory
   - Temporary files may require significant space

## Citation

If you use this tool in your research, please cite:
```
[Citation information to be added]
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please read our contributing guidelines before submitting pull requests.

## Support

For bug reports and feature requests, please use the GitHub issue tracker.

For questions and discussions:
- Open an issue on GitHub
- Contact: [contact information]