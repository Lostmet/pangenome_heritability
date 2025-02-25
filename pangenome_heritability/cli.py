import click
import shutil
import time
import os
from pathlib import Path
from .config import Config
from .variant_processing.vcf_parser import process_variants, filter_vcf
from .variant_processing.fasta_generator import generate_fasta_sequences
from .alignment.mafft_wrapper import run_alignments 
from .kmer.window_generator import process_fasta_files, save_kmer_results_to_csv, process_and_merge_results, parse_fasta_with_metadata
from .kmer.comparison import process_comparison_results
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
    """A Python tool for pangenome heritability analysis."""
    pass
#variants overlap, group and nSV generate
@cli.command("process-vcf")
@click.option('--vcf', required=True, help='Input VCF file')
@click.option('--ref', required=True, help='Reference FASTA file')
@click.option('--out', required=True, help='Output directory for processed variants and FASTA')
@click.option('--threads', default=8, help='Number of threads')
def process_vcf(vcf: str, ref: str, out: str, threads: int):
    """Variants overlap, group and nSV generate"""
    try:
        start_time = time.time()
        click.echo("Processing VCF and generating nSV.vcf...")
        config = Config(vcf_file=vcf, ref_fasta=ref, output_dir=out, threads = threads)
        grouped_variants_list, var_bp_all, var_bp_max, single_sv_count, \
        multi_bp, multi_bp_max, percentage_sv_overlapped, \
        variant_max, single_group= process_variants(config) # var_bp是ref-alt的绝对值之和，不考虑inv，只多不少
        click.echo(f"Total base pairs to be processed: {var_bp_all:,}, max per variant: {var_bp_max:,}")
        click.echo(f'Base pairs to be processed after grouping: {multi_bp:,}, max per variant: {multi_bp_max:,}')
        click.echo(f'Max variant: chromosome = {variant_max.chrom}, position = {variant_max.start:,}')
        click.echo(f'{single_sv_count:,} SVs excluded due to no overlap')
        click.echo(f'Percentage of overlapping SVs: {percentage_sv_overlapped:.2f}%')
        click.echo('Processing SVs no overlapped and...')
        filter_vcf(config, single_group)  
        end_time = time.time()
        total_time = end_time - start_time
        hours = total_time // 3600  
        minutes = (total_time % 3600) // 60 
        seconds = total_time % 60  
        print(f"Total time taken: {int(hours)}:{int(minutes):02d}:{int(seconds):02d}")
    except Exception as e:
        click.echo(f"Error in process-vcf: {str(e)}", err=True)
        raise click.Abort()


# one line
@cli.command("run-all")
@click.option('--vcf', required=True, help='Input VCF file')
@click.option('--ref', required=True, help='Reference FASTA file')
@click.option('--out', required=True, help='Output directory')
@click.option('--threads', default=10, type=int, help='Number of threads')

def run_all(vcf: str, ref: str, out: str, threads: int):
    """Run the complete pipeline."""
    try:
        start_time = time.time()
        click.echo("Step 1: Processing VCF and generating FASTA...")
        config = Config(vcf_file=vcf, ref_fasta=ref, output_dir=out, threads = threads)
        grouped_variants_list, var_bp_all, var_bp_max, single_sv_count, \
        multi_bp, multi_bp_max, percentage_sv_overlapped, \
        variant_max, single_group= process_variants(config) # var_bp是ref-alt的绝对值之和，不考虑inv，只多不少

        click.echo(f"Total base pairs to be processed: {var_bp_all:,}, max per variant: {var_bp_max:,}")
        # 用的最多的应该是实际用到的bp：multi_bp
        click.echo(f'Base pairs to be processed after grouping: {multi_bp:,}, max per variant: {multi_bp_max:,}')
        # 输出一下variant_max
        click.echo(f'Max variant: chromosome = {variant_max.chrom}, position = {variant_max.start:,}')
        # 输出一下有多少SV是没重叠的
        click.echo(f'{single_sv_count:,} SVs excluded due to no overlap')
        # 输出重叠的SV占所有的SV的比例
        click.echo(f'Percentage of overlapping SVs: {percentage_sv_overlapped:.2f}%')
        # 粗略估算预计用时
        # 将总用时转换为小时、分钟、秒
        est_time = multi_bp / 3000
        hours = est_time // 3600  # 总小时数
        minutes = (est_time % 3600) // 60  # 总分钟数
        # 格式化输出
        click.echo(f"Estimated time required: {int(hours)}:{int(minutes):02d}:00")

        ## 开始生成非重叠的vcf文件
            # 输出未重叠的vcf文件
        click.echo('Processing SVs no overlapped and...')
        filter_vcf(config, single_group)  

        grouped_variants_dict = {}
        for group in grouped_variants_list:
            grouped_variants_dict.setdefault(group.chrom, []).append(group) # 得到{"1":[], "2":[]} chrom_number
        total_groups = sum(len(groups) for groups in grouped_variants_dict.values()) 
        click.echo(f"Overlapping variants have been grouped, Group count: {total_groups:,}")
        #print(f'grouped_variants_dict:{grouped_variants_dict}')

        # ✅ 解包生成的输出，获取 fasta 文件路径和 has_insertion
        fasta_path, has_insertion_dict, poly_ins_list = generate_fasta_sequences(config, grouped_variants_dict, total_groups)
        #print(f"打印一下list：{poly_ins_list}，还有dict：{has_insertion_dict}")
        click.echo(f"FASTA file generated: {fasta_path}")


        click.echo("Step 2: Running alignments...")
        alignments_config = Config(
            grouped_variants_file=fasta_path,  # 只传递路径
            ref_fasta=ref,
            output_dir=out,
            threads=threads
        )

        # ✅ 传递参数，路径和has_insertion_dict,poly_ins_list以判断是否存在poly_ins
        run_alignments(alignments_config, fasta_path, has_insertion_dict, poly_ins_list)
        #print(f"dict:{has_insertion_dict}, list:{poly_ins_list}")
        click.echo(f"Alignments completed. Results saved in alignment_results directory")

        click.echo("Step 3: Processing K-mers...")
        alignments_dir = os.path.join(out, "alignment_results")
        intermediate_csv = os.path.join(out, "comparison_results.csv")
        final_csv = os.path.join(out, "output_final_results.csv")
        
        # 获取genome metadata
        genome_metadata = parse_fasta_with_metadata(fasta_path)
        # print(f"meta: {genome_metadata}")
        click.echo("a. Generating kmers...")
        results = process_fasta_files(
            alignments_dir, 
            has_insertion_dict, 
         #传入has_insertion_dict:{'Group_2_1_906668': True, 'Group_2_2_906670': True, 'Group_2_3_906736': False}
            genome_metadata=genome_metadata,
            max_workers=threads
        )
        click.echo("Kmers generated successfully")
        click.echo("b. Saving kmers to .csv for next processing...")
        save_kmer_results_to_csv(results, intermediate_csv)
        click.echo("Done")
        #，merge步骤
        process_and_merge_results(intermediate_csv, final_csv, threads) 

        click.echo(f"K-mer processing completed. Results saved in {final_csv}")

        click.echo("Step 4: Converting to VCF format and Generating matrix...")
        # 创建用于储存矩阵的文件夹
        matrix_dir = Path(out) / "matrix_results"
        matrix_dir.mkdir(parents=True, exist_ok=True)
        click.echo("a. Generating D matrix...")
        # Step 1: 处理 diff_array，生成 D_matrix
        process_diff_array(final_csv, matrix_dir)
        # Step 2: 读取 VCF，生成 X_matrix，并且储存sample_names
        click.echo("b. Generating X matrix...")
        sample_names = process_vcf_to_x_matrix(vcf, matrix_dir)
        # Step 3: 计算 T_matrix
        click.echo("c. Generating T matrix...")
        compute_t_matrix(matrix_dir)

        click.echo("D, X, T matrix generated")
        click.echo("Generating rSV_meta.csv for meta information of rSVs...(processing pos, ref, alt)")
        save_rSV_meta(final_csv, out, threads)#生成rSV_meta，用于查阅ref，alt等信息
        rSV_meta_csv = os.path.join(out, "rSV_meta.csv")
        gt_matrix = os.path.join(out, "GT_matrix.csv")
        output_vcf = os.path.join(out, "output.vcf")
        rSV_vcf = os.path.join(out, "pangenome_rSV.vcf")
        click.echo("4. Generating GT matrix and VCF files...")
        extract_vcf_sample(rSV_meta_csv, gt_matrix, matrix_dir, threads)#生成gt_matrix，用于填充vcf的GT
        #扔进去sample_names用于生成vcf
        vcf_generate(sample_names, rSV_meta_csv, output_vcf, gt_matrix, rSV_vcf)#最终生成rSV.vcf

        click.echo("pangenome_rSV.vcf generated!")

        #第五步，检测不合理的GT数据
        click.echo("Step 5: Detect abnormal GT data...")
        detect_abnormal(matrix_dir)

        click.echo("All steps completed successfully!")
        end_time = time.time()
        # 计算总用时
        total_time = end_time - start_time
        
        # 将总用时转换为小时、分钟、秒
        hours = total_time // 3600  # 总小时数
        minutes = (total_time % 3600) // 60  # 总分钟数
        seconds = total_time % 60  # 剩余秒数
        
        # 格式化输出
        click.echo(f"Total time taken: {int(hours)}:{int(minutes):02d}:{int(seconds):02d}")
    except RuntimeError as e:
        click.echo(f"Error: {str(e)}", err=True)
        raise click.Abort()
    except Exception as e:
        click.echo(f"Error in run-all: {str(e)}", err=True)
        raise click.Abort()

# 有out_final_resuls.csv的可以用这个接着做
@cli.command("make-meta")
@click.option('--vcf', required=True, help='Input VCF file')
@click.option('--out', required=True, help='Output directory')
@click.option('--threads', default=10, type=int, help='Number of threads')

def make_meta(vcf: str, out: str, threads: int):
    final_csv = os.path.join(out, "output_final_results.csv")
    click.echo("Step 4: Converting to VCF format and Generating matrix...")
    # 用于储存矩阵的文件夹
    matrix_dir = Path(out) / "matrix_results"

    click.echo("a. Generating D matrix...")
    # Step 1: 处理 diff_array，生成 D_matrix
    process_diff_array(final_csv, matrix_dir)
    # Step 2: 读取 VCF，生成 X_matrix，并且储存sample_names
    click.echo("b. Generating X matrix...")
    sample_names = process_vcf_to_x_matrix(vcf, matrix_dir)
    # Step 3: 计算 T_matrix
    click.echo("c. Generating T matrix...")
    compute_t_matrix(matrix_dir)

    click.echo("D, X, T matrix generated")

    click.echo("Generating rSV_meta.csv for meta information of rSVs...(processing pos, ref, alt)")
    save_rSV_meta(final_csv, out, threads)#生成rSV_meta，用于查阅ref，alt等信息
    rSV_meta_csv = os.path.join(out, "rSV_meta.csv")
    gt_matrix = os.path.join(out, "GT_matrix.csv")
    output_vcf = os.path.join(out, "output.vcf")
    rSV_vcf = os.path.join(out, "pangenome_rSV.vcf")
    click.echo("4. Generating GT matrix and VCF files...")
    extract_vcf_sample(rSV_meta_csv, gt_matrix, matrix_dir, threads)#生成gt_matrix，用于填充vcf的GT
    #扔进去sample_names用于生成vcf
    vcf_generate(sample_names, rSV_meta_csv, output_vcf, gt_matrix, rSV_vcf)#最终生成rSV.vcf

    click.echo("pangenome_rSV.vcf generated!")

    #第五步，检测不合理的GT数据
    click.echo("Step 5: Detect abnormal GT data...")
    detect_abnormal(matrix_dir)

    click.echo("All steps completed successfully!")

# 生成了meta之后用这个继续做
@cli.command("make-vcf")
@click.option('--vcf', required=True, help='Input VCF file')
@click.option('--out', required=True, help='Output directory')
@click.option('--threads', default=10, type=int, help='Number of threads')
def make_vcf(vcf, out ,threads):
    sample_names = sample_name_contract(vcf)
    matrix_dir = Path(out) / "matrix_results"
    rSV_meta_csv = os.path.join(out, "rSV_meta.csv")
    gt_matrix = os.path.join(out, "GT_matrix.csv")
    output_vcf = os.path.join(out, "output.vcf")
    rSV_vcf = os.path.join(out, "pangenome_rSV.vcf")
    click.echo("4. Generating GT matrix and VCF files...")
    extract_vcf_sample(rSV_meta_csv, gt_matrix, matrix_dir, threads)#生成gt_matrix，用于填充vcf的GT
    #扔进去sample_names用于生成vcf
    vcf_generate(sample_names, rSV_meta_csv, output_vcf, gt_matrix, rSV_vcf)#最终生成rSV.vcf

    click.echo("pangenome_rSV.vcf generated!")

    #第五步，检测不合理的GT数据
    click.echo("Step 5: Detect abnormal GT data...")
    detect_abnormal(matrix_dir)

    click.echo("All steps completed successfully!")
