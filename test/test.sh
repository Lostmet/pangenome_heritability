panherit process-vcf \
    --vcf path/to/pangenome_heritability/test/test.vcf.gz \
    --ref path/to/pangenome_heritability/test/test.fasta \
    --out path/to/pangenome_heritability/test

panherit run-alignments \
    --grouped-variants path/to/pangenome_heritability/test/variants.fasta \
    --ref path/to/pangenome_heritability/test/test.fasta \
    --out path/to/pangenome_heritability/test \
    --threads 1

panherit process-kmers \
    --alignments path/to/pangenome_heritability/test/alignment_results \
    --window-size 4 \
    --out path/to/pangenome_heritability/test/kmers_directory \
    --threads 1

panherit convert-to-plink \
    --csv-file path/to/pangenome_heritability/test/kmers_directory/output_final_results.csv \
    --grouped-variants path/to/pangenome_heritability/test/variants.fasta \
    --vcf-file path/to/pangenome_heritability/test/test.vcf.gz \
    --output-dir path/to/pangenome_heritability/test/ 



panherit run-all \
    --vcf /storage/yangjianLab/yuanpeixiong/test/pangenome_heritability/test/test.vcf.gz \
    --ref /storage/yangjianLab/yuanpeixiong/test/pangenome_heritability/test/test.fasta \
    --out /storage/yangjianLab/yuanpeixiong/test/pangenome_heritability/test \
    --window-size 4 \
    --threads 4


