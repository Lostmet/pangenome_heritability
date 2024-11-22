import click
from .config import Config
from .variant_processing.vcf_parser import process_variants
from .variant_processing.fasta_generator import generate_fasta_sequences
from .alignment.muscle_wrapper import run_alignments
from .kmer.window_generator import process_windows
from .genotype.plink_converter import convert_to_plink

@click.group()
def cli():
    """A CLI tool for pangenome heritability analysis."""
    pass
#variant overlap and group
@cli.command("process-vcf")
@click.option('--vcf', required=True, help='Input VCF file')
@click.option('--ref', required=True, help='Reference FASTA file')
@click.option('--out', required=True, help='Output directory for processed variants and FASTA')
def process_vcf(vcf: str, ref: str, out: str):
    """Process VCF file, group overlapping variants, and generate FASTA."""
    try:
        # Step 1
        config = Config(vcf_file=vcf, ref_fasta=ref, output_dir=out)
        
        # Step 2
        grouped_variants = process_variants(config)
        click.echo(f"Variants processed and grouped. Total groups: {len(grouped_variants)}")
        
        # Step 3
        fasta_path = generate_fasta_sequences(config, grouped_variants)
        click.echo(f"FASTA file generated: {fasta_path}")
        
    except Exception as e:
        click.echo(f"Error in process-vcf: {str(e)}", err=True)
        raise click.Abort()

# Step 2: Run alignments
@cli.command("run-alignments")
@click.option('--grouped-variants', required=True, help='Grouped variants file')
@click.option('--ref', required=True, help='Reference FASTA file')
@click.option('--out', required=True, help='Output directory for alignments')
@click.option('--threads', default=1, help='Number of threads')
def run_alignments_cmd(grouped_variants: str, ref: str, out: str, threads: int):
    """Run alignments for grouped variants."""
    try:
        config = Config(grouped_variants_file=grouped_variants, ref_fasta=ref, output_dir=out, threads=threads)
        alignments = run_alignments(config, grouped_variants)
        click.echo(f"Alignments completed. Output saved to {out}")
    except Exception as e:
        click.echo(f"Error in run-alignments: {str(e)}", err=True)
        raise click.Abort()

# Step 3: Process K-mer windows
@cli.command("process-kmers")
@click.option('--alignments', required=True, help='Input alignments directory')
@click.option('--window-size', default=4, help='K-mer window size')
@click.option('--out', required=True, help='Output directory for K-mer results')
def process_kmers(alignments: str, window_size: int, out: str):
    """Process K-mer windows based on alignments."""
    try:
        config = Config(alignments_dir=alignments, window_size=window_size, output_dir=out)
        kmer_results = process_windows(config, alignments)
        click.echo(f"K-mers processed. Results saved to {out}")
    except Exception as e:
        click.echo(f"Error in process-kmers: {str(e)}", err=True)
        raise click.Abort()

# Step 4: Convert to PLINK
@cli.command("convert-to-plink")
@click.option('--kmer-results', required=True, help='K-mer results directory')
@click.option('--out', required=True, help='Output directory for PLINK files')
def convert_to_plink_cmd(kmer_results: str, out: str):
    """Convert K-mer results to PLINK binary files."""
    try:
        config = Config(kmer_results_dir=kmer_results, output_dir=out)
        convert_to_plink(config, kmer_results)
        click.echo(f"PLINK files generated. Output saved to {out}")
    except Exception as e:
        click.echo(f"Error in convert-to-plink: {str(e)}", err=True)
        raise click.Abort()