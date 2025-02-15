import click
import shutil
import sys
import os
import pandas as pd
from .config import Config
from .variant_processing.vcf_parser import process_variants
from .variant_processing.fasta_generator import generate_fasta_sequences
from .alignment.muscle_wrapper import run_alignments
from .kmer.window_generator import process_fasta_files, save_kmer_results_to_csv, process_and_merge_results, read_fasta_files , parse_fasta_with_metadata
from .kmer.comparison import process_comparison_results
from .genotype.genotype_mapper import create_ped_and_map_files, create_vcf_file, process_diff_array, process_vcf_to_x_matrix, compute_t_matrix, save_rsv_meta, extract_vcf_sample, vcf_generate, detect_abnormal
from .utils.logging_utils import get_logger
from .variant_processing.vcf_parser import VariantGroup

def check_tools(*tools):
    """
    Check if required tools are available in the user's PATH.
    Raise an error if any tool is missing.
    """
    missing_tools = [tool for tool in tools if shutil.which(tool) is None]
    if missing_tools:
        raise RuntimeError(
            f"The following tools are missing from your environment: {', '.join(missing_tools)}. "
            f"Please ensure they are installed and accessible from your PATH."
        )

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

        # ✅ 先调用 generate_fasta_sequences()，正确解包返回的 Tuple
        fasta_file, has_insertion = generate_fasta_sequences(config, grouped_variants)
        print(f"DEBUG: fasta_file = {fasta_file}")  # 确保是文件路径字符串
        print(f"DEBUG: Is fasta_file a string? {isinstance(fasta_file, str)}")


        # ✅ 传入两个参数，防止参数不匹配
        alignments = run_alignments(config, fasta_file, has_insertion)

        click.echo(f"Alignments completed. Output saved to {out}")
    except RuntimeError as e:
        click.echo(f"Error: {str(e)}", err=True)
        raise click.Abort()
    except Exception as e:
        click.echo(f"Error in run-alignments: {str(e)}", err=True)
        raise click.Abort()


@cli.command("process-kmers")
@click.option('--alignments', required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True), 
              help='Input alignments directory containing FASTA files')
@click.option('--window-size', default=4, type=int, help='K-mer window size (default: 4)')
@click.option('--grouped-variants', required=True, type=click.Path(exists=True, file_okay=True, dir_okay=False), 
              help='Path to the FASTA file containing grouped variants.')
@click.option('--out', required=True, type=click.Path(file_okay=False, dir_okay=True), 
              help='Output directory for final results')
@click.option('--threads', type=int, default=None, help='Maximum number of worker processes (default: CPU count - 1)')
def process_kmers(alignments: str, window_size: int, grouped_variants: str, out: str, threads: int):
    """
    Process K-mer windows, including:
    1. Parse genome FASTA metadata
    2. Process alignment results to generate K-mer comparison results
    3. Merge adjacent k-mer windows and output results
    """
    try:
        # Ensure output directory exists
        os.makedirs(out, exist_ok=True)

        # Define output file paths
        intermediate_csv = os.path.join(out, "comparison_results.csv")
        final_csv = os.path.join(out, "output_final_results.csv")

        # Step 1: Parse genome FASTA metadata
        click.echo(f"Step 1: Parsing genome FASTA metadata: {grouped_variants}")
        genome_metadata = parse_fasta_with_metadata(grouped_variants)
        click.echo("Genome FASTA metadata parsing completed")

        # Step 2: Process alignment results
        click.echo(f"Step 2: Processing FASTA files with window size {window_size}")
        results = process_fasta_files(
            alignments, 
            genome_metadata=genome_metadata, 
            k=window_size, 
            max_workers=threads
        )
        save_kmer_results_to_csv(results, intermediate_csv)
        click.echo(f"K-mer comparison results saved to {intermediate_csv}")

        # Step 3: Process and merge results
        click.echo("Step 3: Merging adjacent k-mer windows")
        process_and_merge_results(intermediate_csv, final_csv)
        click.echo(f"Final results saved to {final_csv}")

        click.echo("K-mer processing completed!")
    except Exception as e:
        click.echo(f"Error during processing: {str(e)}", err=True)
        raise click.Abort()


# Step 4: Convert to PLINK
@cli.command("convert-to-plink")
@click.option('--csv-file', required=True, help='Path to CSV file containing final comparison results and metadata')
@click.option('--vcf-file', required=True, help='Path to original VCF file')
@click.option('--output-dir', required=True, help='Output directory for PLINK files')
def convert_to_plink_cmd(csv_file: str, vcf_file: str, output_dir: str):
    """
    Generate PLINK files using comparison results and metadata, marking SVs and RSVs.
    """
    try:
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate output file prefix
        output_prefix = os.path.join(output_dir, "pangenome")
        
        # Read CSV file
        csv_data = pd.read_csv(csv_file)
        
        # Modified to add vcf_file parameter
        create_ped_and_map_files(
            csv_data=csv_data,
            vcf_file=vcf_file,  # Add this parameter
            output_prefix=output_prefix
        )
        
        click.echo(f"PLINK files successfully generated in {output_dir}")
        click.echo(f"- PED file: {output_prefix}.ped")
        click.echo(f"- MAP file: {output_prefix}.map")

    except Exception as e:
        click.echo(f"Error generating PLINK files: {str(e)}", err=True)
        raise click.Abort()

# one line
@cli.command("run-all")
@click.option('--vcf', required=True, help='Input VCF file')
@click.option('--ref', required=True, help='Reference FASTA file')
@click.option('--out', required=True, help='Output directory')
@click.option('--window-size', default=4, type=int, help='K-mer window size (default: 4)')
@click.option('--threads', default=1, type=int, help='Number of threads')
def run_all(vcf: str, ref: str, out: str, window_size: int, threads: int):
    """Run the complete pipeline."""
    try:
        click.echo("Step 1: Processing VCF and generating FASTA...")
        config = Config(vcf_file=vcf, ref_fasta=ref, output_dir=out)
        grouped_variants_list = process_variants(config)
        grouped_variants_dict = {}
        for group in grouped_variants_list:
            grouped_variants_dict.setdefault(group.chrom, []).append(group)
        
        # ✅ 解包生成的输出，获取 fasta 文件路径和 has_insertion
        fasta_path, has_insertion_dict, poly_ins_list = generate_fasta_sequences(config, grouped_variants_dict)

        click.echo(f"FASTA file generated: {fasta_path}")

        click.echo("Step 2: Running alignments...")
        alignments_config = Config(
            grouped_variants_file=fasta_path,  # 只传递路径
            ref_fasta=ref,
            output_dir=out,
            threads=threads
        )

        # ✅ 传递两个参数，路径和has_insertion
        run_alignments(alignments_config, fasta_path, has_insertion_dict, poly_ins_list)
        click.echo(f"Alignments completed. Results saved in alignment_results directory")

        click.echo("Step 3: Processing K-mers...")
        alignments_dir = os.path.join(out, "alignment_results")
        intermediate_csv = os.path.join(out, "comparison_results.csv")
        processed_csv = os.path.join(out, "processed_comparison_results.csv")
        final_csv = os.path.join(out, "output_final_results.csv")
        
        # 获取genome metadata
        genome_metadata = parse_fasta_with_metadata(fasta_path)
        
        results = process_fasta_files(
            alignments_dir, 
            genome_metadata=genome_metadata,
            k=window_size, 
            max_workers=threads
        )
        save_kmer_results_to_csv(results, intermediate_csv)
        process_comparison_results(intermediate_csv, processed_csv)
        process_and_merge_results(processed_csv, final_csv)
        click.echo(f"K-mer processing completed. Results saved in {final_csv}")

        click.echo("Step 4: Converting to VCF format...")
        vcf_prefix = os.path.join(out, "pangenome")
        #csv_data = pd.read_csv(final_csv)
        # Step 1: 处理 diff_array，生成 D_matrix
        process_diff_array(final_csv, out)
        # Step 2: 读取 VCF，生成 X_matrix
        process_vcf_to_x_matrix(vcf, out)
        # Step 3: 计算 T_matrix
        compute_t_matrix(out)

        click.echo("D, X, T matrix generated")

        save_rsv_meta(final_csv, out)#生成rsv_meta，用于查阅ref，alt等信息
        rsv_meta_csv = os.path.join(out, "rsv_meta.csv")
        gt_matrix = os.path.join(out, "GT_matrix.csv")
        output_vcf = os.path.join(out, "output.vcf")
        rsv_vcf = os.path.join(out, "pangenome_rsv.vcf")
        extract_vcf_sample(rsv_meta_csv, gt_matrix, out)#生成gt_matrix，用于填充vcf的GT
        vcf_generate(vcf, rsv_meta_csv, output_vcf, gt_matrix, rsv_vcf)#最终生成rsv.vcf


        click.echo("pangenome_rsv.vcf generated!")

        #第五步，检测不合理的GT数据
        click.echo("Step 5: Detect abnormal GT data...")
        detect_abnormal(out)

        click.echo("All steps completed successfully!")
    except RuntimeError as e:
        click.echo(f"Error: {str(e)}", err=True)
        raise click.Abort()
    except Exception as e:
        click.echo(f"Error in run-all: {str(e)}", err=True)
        raise click.Abort()


@cli.command("convert-to-vcf")
@click.option('--csv-file', required=True, help='Path to CSV file containing final comparison results and metadata')
@click.option('--vcf-file', required=True, help='Path to original VCF file')
@click.option('--grouped-variants', required=True, help='Path to FASTA file containing variant information')
@click.option('--output-dir', required=True, help='Output directory')
def convert_to_vcf_cmd(csv_file: str, vcf_file: str, grouped_variants: str, output_dir: str):
    """Generate VCF file using comparison results and metadata."""
    try:
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate output file prefix
        output_prefix = os.path.join(output_dir, "pangenome")
        
        # Read CSV file
        csv_data = pd.read_csv(csv_file)
        
        # Generate VCF file
        create_vcf_file(
            csv_data=csv_data,
            vcf_file=vcf_file,
            grouped_variants=grouped_variants,  # Add new parameter
            output_prefix=output_prefix
        )
        
        click.echo(f"VCF file successfully generated in {output_dir}")

    except Exception as e:
        click.echo(f"Error generating VCF file: {str(e)}", err=True)
        raise click.Abort()
