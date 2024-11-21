import click
from .config import Config
from .variant_processing.vcf_parser import process_variants
from .alignment.muscle_wrapper import run_alignments
from .kmer.window_generator import process_windows
from .genotype.plink_converter import convert_to_plink

@click.command()
@click.option('--vcf', required=True, help='Input VCF file')
@click.option('--ref', required=True, help='Reference FASTA file')
@click.option('--out', required=True, help='Output directory')
@click.option('--threads', default=1, help='Number of threads')
@click.option('--window-size', default=4, help='K-mer window size')
def main(vcf: str, ref: str, out: str, threads: int, window_size: int):
    """Generate PLINK binary files from structural variants"""
    try:
        config = Config(
            vcf_file=vcf,
            ref_fasta=ref,
            output_dir=out,
            threads=threads,
            window_size=window_size
        )
        
        grouped_variants = process_variants(config)
        alignments = run_alignments(config, grouped_variants)
        kmer_results = process_windows(config, alignments)
        convert_to_plink(config, kmer_results)
        
    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        raise click.Abort()