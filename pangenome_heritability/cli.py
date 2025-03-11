import click
import shutil
import time
import os
import logging
from pathlib import Path
from .config import Config
from .variant_processing.vcf_parser import process_variants, filter_vcf
from .variant_processing.fasta_generator import generate_fasta_sequences
from .alignment.mafft_wrapper import run_alignments 
from .rSV.window_generator import process_fasta_files, save_rSV_results_to_csv, process_and_merge_results, parse_fasta_with_metadata, nrSV_vcf_generate
from .rSV.comparison import process_comparison_results
from .genotype.genotype_mapper import process_diff_array, process_vcf_to_x_matrix, compute_t_matrix, save_rSV_meta, extract_vcf_sample, vcf_generate, detect_abnormal, sample_name_contract
from .utils.logging_utils import get_logger

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
    """A Python tool for refined Structured Variants analysis."""
    pass

@cli.command("process-vcf")
@click.option('--vcf', required=True, help='Input VCF file')
@click.option('--ref', required=True, help='Reference FASTA file')
@click.option('--cutoff', default=0.9, type=float, help="Threshold for assigning poly-alts to separate rSVs.")
@click.option('--out', required=True, help='Output directory for processed variants and FASTA files')
@click.option('--threads', default=8, type=int, help='Number of threads')

def process_vcf(vcf: str, ref: str, cutoff, out: str, threads: int):
    """Process overlapping variants, group them, and generate nSVs. 
    This step produces 'variants_pre_aligned.fasta' and 'nSV.vcf'."""
    try:
        start_time = time.time()
        click.echo("[Step 1] Processing VCF and generating FASTA...")
        config = Config(vcf_file=vcf, ref_fasta=ref, output_dir=out, threads=threads)
        _, var_bp_all, var_bp_max, single_sv_count, \
        multi_bp, multi_bp_max, percentage_sv_overlapped, \
        variant_max, single_group, _, _ = process_variants(config)

        click.echo(f"Total base pairs to be processed: {var_bp_all:,}, max per variant: {var_bp_max:,}")
        click.echo(f"After grouping: {multi_bp:,} base pairs, max per variant: {multi_bp_max:,}")
        click.echo(f"Max variant located at chromosome {variant_max.chrom}, position {variant_max.start:,}")
        click.echo(f"{single_sv_count:,} SVs excluded due to lack of overlap")
        click.echo(f"Percentage of overlapping SVs: {percentage_sv_overlapped:.2f}%")
        
        click.echo("Processing non-overlapping SVs (nSVs) and generating nSV.vcf...")
        _ = filter_vcf(config, single_group)
        
        end_time = time.time()
        total_time = end_time - start_time
        hours, remainder = divmod(total_time, 3600)
        minutes, seconds = divmod(remainder, 60)
        
        logging.basicConfig(
            filename=os.path.join(out, f"{out}_process-vcf.log"),
            level=logging.INFO,
            format="%(asctime)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        log_messages = [
            f"Total runtime: {int(hours)}:{int(minutes):02d}:{int(seconds):02d}",
            f"Total base pairs processed: {var_bp_all:,}",
            f"Max base pairs per variant: {var_bp_max:,}",
            f"Grouped base pairs: {multi_bp:,}, max per variant: {multi_bp_max:,}",
            f"Max variant: chromosome {variant_max.chrom}, position {variant_max.start:,}",
            f"Excluded SVs (no overlap): {single_sv_count:,}",
            f"Overlapping SVs percentage: {percentage_sv_overlapped:.2f}%"
        ]
        
        for msg in log_messages:
            logging.info(msg)
    except Exception as e:
        click.echo(f"Error in process-vcf: {str(e)}", err=True)
        logging.error(f"Error in process-vcf: {str(e)}")
        raise click.Abort()



@cli.command("run-all")
@click.option('--vcf', required=True, help='Input VCF file')
@click.option('--ref', required=True, help='Reference FASTA file')
@click.option('--cutoff', default=0.9, type=float, help="Threshold for separating poly-alts into distinct rSVs.")
@click.option('--out', required=True, help='Output directory')
@click.option('--threads', default=10, type=int, help='Number of threads')

def run_all(vcf: str, ref: str, cutoff: float, out: str, threads: int):
    """Execute the complete variant processing pipeline."""
    try:
        start_time = time.time()
        click.echo("[Step 1] Processing VCF and generating FASTA...")
        config = Config(vcf_file=vcf, ref_fasta=ref, output_dir=out, threads=threads)
        grouped_variants_list, var_bp_all, var_bp_max, single_sv_count, \
        multi_bp, multi_bp_max, percentage_sv_overlapped, \
        variant_max, single_group, inv_count, variant_count = process_variants(config)

        click.echo(f"Total base pairs to be processed: {var_bp_all:,}, max per variant: {var_bp_max:,}")
        click.echo(f"After grouping: {multi_bp:,} base pairs, max per variant: {multi_bp_max:,}")
        click.echo(f"Max variant located at: chromosome {variant_max.chrom}, position {variant_max.start:,}")
        click.echo("If the runtime is too long, you can remove this variant to speed up the process.")
        click.echo(f"{single_sv_count:,} SVs excluded due to lack of overlap")
        click.echo(f"Percentage of overlapping SVs: {percentage_sv_overlapped:.2f}%")
        
        est_time = multi_bp / 3000
        hours, remainder = divmod(est_time, 3600)  
        minutes, _ = divmod(remainder, 60)  
        click.echo(f"Estimated runtime: {int(hours)}:{int(minutes):02d}:00")
        click.echo(click.style("Note: Actual runtime depends on CPU performance.", fg="yellow"))

        click.echo("Processing non-overlapping SVs (nSVs) and generating nSV.vcf...")
        nSV_name = filter_vcf(config, single_group)
        
        grouped_variants_dict = {}
        for group in grouped_variants_list:
            grouped_variants_dict.setdefault(group.chrom, []).append(group)
        total_groups = sum(len(groups) for groups in grouped_variants_dict.values())
        click.echo(f"Overlapping variants grouped: {total_groups:,} groups")

        fasta_path, has_insertion_dict, poly_ins_list = generate_fasta_sequences(config, grouped_variants_dict, total_groups)
        click.echo(f"FASTA file created: {fasta_path}")

        click.echo("[Step 2] Running alignments...")
        alignments_config = Config(grouped_variants_file=fasta_path, ref_fasta=ref, output_dir=out, threads=threads)
        run_alignments(alignments_config, fasta_path, has_insertion_dict, poly_ins_list)
        click.echo("Alignments completed. Results saved in 'alignment_results' directory.")

        click.echo("[Step 3] Generating rSVs...")
        alignments_dir = os.path.join(out, "alignment_results")
        intermediate_csv = os.path.join(out, "Scanning_results.csv")
        final_csv = os.path.join(out, "merged_rSVs.csv")
        nrSV_csv = os.path.join(out, "nrSV_meta.csv")

        genome_metadata = parse_fasta_with_metadata(fasta_path)
        click.echo("Scanning aligned FASTA files...")
        results = process_fasta_files(alignments_dir, has_insertion_dict, genome_metadata, max_workers=threads)
        click.echo("All FASTA files successfully scanned.")
        
        click.echo(f"Temporarily caching results in {intermediate_csv} for further processing...")
        save_rSV_results_to_csv(results, intermediate_csv)
        click.echo("Done.")

        nrSV_list = process_and_merge_results(intermediate_csv, final_csv, threads, has_insertion_dict, cutoff, nrSV_csv)
        click.echo(f"rSV processing completed. Results saved in {final_csv}")
        os.remove(intermediate_csv)
        
        click.echo("[Step 4] Converting to VCF format and generating matrices...")
        if nrSV_list:
            nrSV_count = nrSV_vcf_generate(nrSV_csv, config, nrSV_list, nSV_name)
        
        click.echo("Generating rSV_meta.csv with meta information (positions, reference, alternate alleles)...")
        save_rSV_meta(final_csv, out, threads)

        matrix_dir = Path(out) / "matrix_results"
        matrix_dir.mkdir(parents=True, exist_ok=True)
        click.echo("Generating matrices...")
        click.echo("Creating D matrix...")
        process_diff_array(final_csv, matrix_dir)
        click.echo("Creating X matrix...")
        sample_names = process_vcf_to_x_matrix(vcf, matrix_dir)
        click.echo("Creating T matrix...")
        compute_t_matrix(matrix_dir)
        click.echo("Matrices successfully generated in 'matrix_results'.")

        rSV_meta_csv = os.path.join(out, "rSV_meta.csv")
        gt_matrix = os.path.join(out, "GT_matrix.csv")
        output_vcf = os.path.join(out, "rSV_without_genotype.vcf")
        rSV_vcf = os.path.join(out, "rSV.vcf")

        click.echo("Generating GT matrix and rSV.vcf...")
        rSV_count, total_groups = extract_vcf_sample(rSV_meta_csv, gt_matrix, matrix_dir, threads)
        click.echo(f"rSV count: {rSV_count:,}")
        vcf_generate(sample_names, rSV_meta_csv, output_vcf, gt_matrix, rSV_vcf)

        try:
            os.remove(output_vcf)
            os.remove(gt_matrix)
        except:
            pass
        click.echo(f"{rSV_vcf} successfully generated!")
        click.echo("Pipeline execution completed!")

        end_time = time.time()
        total_time = end_time - start_time
        hours, minutes, seconds = total_time // 3600, (total_time % 3600) // 60, total_time % 60

        logging.basicConfig(
            filename=os.path.join(out, f"{out}_run-all.log"),
            level=logging.INFO,
            format="%(asctime)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

        log_messages = [
            f"Total runtime: {int(hours)}:{int(minutes):02d}:{int(seconds):02d}",
            f"Total variants: {variant_count:,}",
            f"INV count: {inv_count:,}",
            f"nSV count: {single_sv_count:,}",
            f"Overlapping SVs: {(variant_count - single_sv_count):,}",
            f"Overlap percentage: {percentage_sv_overlapped:.2f}%",
            f"Total variant groups: {total_groups:,}",
            f"Final rSV count: {rSV_count:,}"
        ]

        for msg in log_messages:
            click.echo(msg)
            logging.info(msg)
    
    except RuntimeError as e:
        click.echo(f"Error: {str(e)}", err=True)
        raise click.Abort()
    except Exception as e:
        click.echo(f"Critical error in run-all: {str(e)}", err=True)
        raise click.Abort()

@cli.command("make-meta")
@click.option('--vcf', required=True, help='Input VCF file')
@click.option('--out', required=True, help='Output directory')
@click.option('--ref', default=None, help='Reference FASTA file')
@click.option('--cutoff', default=0.9, type=float, help="Threshold for assigning poly-alts to separate rSVs.")
@click.option('--threads', default=10, type=int, help='Number of threads')

def make_meta(vcf: str, ref, _, out: str, threads: int):
    '''If you already have merged_rSVs.csv and the 'run-all' step was interrupted at this stage, you can resume from here. 
    Ensure the output path remains unchanged.'''
    start_time = time.time()
    click.echo("[Step 1] Extracting SV information...")
    config = Config(vcf_file=vcf, ref_fasta=ref, output_dir=out, threads=threads)
    _, _, _, single_sv_count, \
    _, _, percentage_sv_overlapped, \
    _, _, inv_count, variant_count = process_variants(config)

    final_csv = os.path.join(out, "merged_rSVs.csv")
    click.echo("[Step 2] Converting to VCF format and generating matrices...")
    
    # Create a directory for storing matrices
    matrix_dir = Path(out) / "matrix_results"
    matrix_dir.mkdir(parents=True, exist_ok=True)
    
    sample_names = sample_name_contract(vcf)
    
    click.echo("1. Generating rSV_meta.csv with metadata (positions, reference, and alternate alleles)...")
    save_rSV_meta(final_csv, out, threads)  # Generates rSV_meta for reviewing ref, alt, etc.

    click.echo("2. Creating matrices...")
    click.echo("a. Generating D matrix...")
    process_diff_array(final_csv, matrix_dir)
    
    click.echo("b. Generating X matrix...")
    sample_names = process_vcf_to_x_matrix(vcf, matrix_dir)
    
    click.echo("c. Generating T matrix...")
    compute_t_matrix(matrix_dir)
    click.echo("D, X, and T matrices successfully generated in 'matrix_results' directory.")
    
    rSV_meta_csv = os.path.join(out, "rSV_meta.csv")
    gt_matrix = os.path.join(out, "GT_matrix.csv")
    output_vcf = os.path.join(out, "rSV_without_genotype.vcf")
    rSV_vcf = os.path.join(out, "rSV.vcf")
    
    click.echo("3. Generating GT matrix and rSV.vcf...")
    rSV_count, total_groups = extract_vcf_sample(rSV_meta_csv, gt_matrix, matrix_dir, threads)  # Generates GT matrix for VCF
    click.echo(f"Total rSV count: {rSV_count:,}")
    
    vcf_generate(sample_names, rSV_meta_csv, output_vcf, gt_matrix, rSV_vcf)  # Generates final rSV.vcf
    
    try:
        os.remove(output_vcf)
        os.remove(gt_matrix)
    except:
        pass
    click.echo(f"{rSV_vcf} successfully generated!")
    click.echo("All steps completed successfully!")
    
    end_time = time.time()
    total_time = end_time - start_time
    hours, remainder = divmod(total_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    logging.basicConfig(
        filename=os.path.join(out, f"{out}_make-meta.log"),
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    log_messages = [
        f"Total runtime: {int(hours)}:{int(minutes):02d}:{int(seconds):02d}",
        f"Total variants: {variant_count:,}",
        f"INV count: {inv_count:,}",
        f"nSV count: {single_sv_count:,}",
        f"Overlapping SVs: {(variant_count - single_sv_count):,}",
        f"Overlap percentage: {percentage_sv_overlapped:.2f}%",
        f"Total variant groups: {total_groups:,}",
        f"Final rSV count: {rSV_count:,}"
    ]
    
    for msg in log_messages:
        click.echo(msg)
        logging.info(msg)
