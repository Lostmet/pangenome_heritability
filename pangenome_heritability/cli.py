import click
import os

from .config import Config
from .variant_processing.vcf_parser import process_variants
from .variant_processing.fasta_generator import generate_fasta_sequences
from .alignment.muscle_wrapper import run_alignments
from .kmer.window_generator import process_fasta_files, save_kmer_results_to_csv, process_chromosome_groups, process_and_merge_results, read_fasta_files
from .kmer.comparison import process_comparison_results
from .genotype.genotype_mapper import convert_to_plink_with_variants
from .utils.logging_utils import get_logger

logger = get_logger(__name__)
@click.group()
def cli():
    """A Python tool for pangenome heritability analysis."""
    pass
#Step1: variant overlap and group
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

#Step2: Run alignments
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
@cli.command("process-kmers")
@click.option('--alignments', required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True), help='Input alignments directory containing FASTA files')
@click.option('--window-size', default=4, type=int, help='K-mer window size (default: 4)')
@click.option('--out', required=True, type=click.Path(file_okay=False, dir_okay=True), help='Output directory for final results')
@click.option('--threads', type=int, default=None, help='Maximum number of worker processes (default: CPU count - 1)')
def process_kmers(alignments: str, window_size: int, out: str, threads: int):
    """
    Process K-mer windows based on alignments, including:
    1. Processing alignments into K-mer comparison results.
    2. Removing redundancy and handling invalid data.
    3. Merging results and removing collinearity.
    """
    try:
        # Ensure output directory exists
        os.makedirs(out, exist_ok=True)

        # Define paths for intermediate and final output files
        intermediate_csv = os.path.join(out, "comparison_results.csv")
        processed_csv = os.path.join(out, "processed_comparison_results.csv")
        final_csv = os.path.join(out, "output_final_results.csv")

        # Step 1: Process alignments into K-mer comparison results
        click.echo(f"Step 1: Processing FASTA files from {alignments} with window size {window_size}")
        results = process_fasta_files(alignments, k=window_size, max_workers=threads)
        save_kmer_results_to_csv(results, intermediate_csv)
        click.echo(f"K-mer comparison results saved to {intermediate_csv}")

        # Step 2: Remove redundancy and handle invalid data
        click.echo("Step 2: Processing chromosome groups and removing redundancy")
        process_comparison_results(intermediate_csv, processed_csv)
        click.echo(f"Processed group results saved to {processed_csv}")

        # Step 3: Merge results and remove collinearity
        click.echo("Step 3: Merging results and removing collinearity")
        process_and_merge_results(processed_csv, final_csv)
        click.echo(f"Final merged results saved to {final_csv}")

        click.echo("K-mer processing complete!")
    except Exception as e:
        click.echo(f"Error in process-kmers: {str(e)}", err=True)
        raise click.Abort()


# Step 4: Convert to PLINK
@cli.command("convert-to-plink")
@click.option('--csv-file', required=True, help='Path to the processed comparison CSV file.')
@click.option('--grouped-variants', required=True, help='Path to the FASTA file containing grouped variants.')
@click.option('--output-dir', required=True, help='Output directory for PLINK files.')
def convert_to_plink_cmd(csv_file: str, fasta_file: str, output_dir: str):
    """CLI command: Replace variant names and generate PLINK files."""
    try:
        # Configure the Config object
        config = Config(
            output_dir=output_dir,
            grouped_variants_file=csv_file,
            ref_fasta=fasta_file
        )

        # Execute main process
        convert_to_plink_with_variants(config)
        click.echo(f"PLINK files successfully generated in {output_dir}")

    except Exception as e:
        click.echo(f"Error in convert-to-plink: {str(e)}", err=True)
        raise click.Abort()

if __name__ == "__main__":
    cli()