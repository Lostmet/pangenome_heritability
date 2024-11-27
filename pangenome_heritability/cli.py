import click
import os

from .config import Config
from .variant_processing.vcf_parser import process_variants
from .variant_processing.fasta_generator import generate_fasta_sequences
from .alignment.muscle_wrapper import run_alignments
from .kmer.window_generator import process_fasta_files, save_kmer_results_to_csv
from .genotype.plink_converter import convert_to_plink

@click.group()
def cli():
    """A Python tool for pangenome heritability analysis."""
    pass
#variant overlap and group
@cli.command("process-vcf")
@click.option('--vcf', required=True, help='Input VCF file')
@click.option('--ref', required=True, help='Reference FASTA file')
@click.option('--out', required=True, help='Output directory for processed variants and FASTA')
def process_vcf(vcf: str, ref: str, out: str):
    """Group overlapping variants, and generate FASTA."""
    try:
        # Step 1: Configuration
        config = Config(vcf_file=vcf, ref_fasta=ref, output_dir=out)
        
        # Step 2: Process variants
        grouped_variants_list = process_variants(config)
        click.echo(f"Variants processed and grouped. Total groups: {len(grouped_variants_list)}")
        
        # Step 2.1: Group variants by chromosome
        grouped_variants_dict = {}
        for group in grouped_variants_list:
            if group.chrom not in grouped_variants_dict:
                grouped_variants_dict[group.chrom] = []
            grouped_variants_dict[group.chrom].append(group)
        
        # Step 3: Generate FASTA file
        fasta_path = generate_fasta_sequences(config, grouped_variants_dict)
        click.echo(f"FASTA file generated: {fasta_path}")
        
    except Exception as e:
        click.echo(f"Error in process-vcf: {str(e)}", err=True)
        raise click.Abort()

# Run alignments
@cli.command("run-alignments")
@click.option('--grouped-variants', required=True, help='Grouped variants file')
@click.option('--ref', required=True, help='Reference FASTA file')
@click.option('--out', required=True, help='Output directory for alignments')
@click.option('--threads', default=1, help='Number of threads')
def run_alignments_cmd(grouped_variants: str, ref: str, out: str, threads: int):
    """Run alignments for grouped variants."""
    try:
        
        config = Config(
            grouped_variants_file=grouped_variants,
            ref_fasta=ref,
            output_dir=out,
            threads=threads
        )
    
        alignments = run_alignments(config, grouped_variants)
        click.echo(f"Alignments completed. Output saved to {out}")
    except Exception as e:
        click.echo(f"Error in run-alignments: {str(e)}", err=True)
        raise click.Abort()

# Step 3: Process K-mer windows
@click.command("process-kmers")
@click.option('--alignments', required=True, help='Input alignments directory')
@click.option('--window-size', default=4, help='K-mer window size')
@click.option('--out', required=True, help='Output directory for K-mer results')
def process_kmers(alignments: str, window_size: int, out: str):
    """Process K-mer windows based on alignments."""
    try:
        if not os.path.exists(alignments):
            raise FileNotFoundError(f"Alignments directory not found: {alignments}")

        logger.info(f"Processing alignments from {alignments} with window size {window_size}")
        results = process_fasta_files(alignments, window_size)

        output_file = os.path.join(out, "kmer_results.csv")
        save_kmer_results_to_csv(results, output_file)

        click.echo(f"K-mers processed. Results saved to {output_file}")
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
if __name__ == "__main__":
    cli()