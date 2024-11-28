panherit process-vcf --vcf /storage/yangjianLab/yuanpeixiong/SV_project/pangenome_heritability/tests/test_data/test.vcf.gz \
                     --ref /storage/yangjianLab/yuanpeixiong/SV_project/pangenome_heritability/tests/test_data/SL5.0.fasta \
                     --out /storage/yangjianLab/yuanpeixiong/pangenome_heritability/test

panherit run-alignments \
    --grouped-variants /storage/yangjianLab/yuanpeixiong/pangenome_heritability/test/variants.fasta \
    --ref /storage/yangjianLab/yuanpeixiong/SV_project/pangenome_heritability/tests/test_data/SL5.0.fasta \
    --out /storage/yangjianLab/yuanpeixiong/SV_project/sv_heritiability/test_data/alignment \
    --threads 20

panherit process-kmers \
    --alignments /storage/yangjianLab/yuanpeixiong/SV_project/sv_heritiability/test_data/alignment/temp_alignments \
    --window-size 4 \
    --out /storage/yangjianLab/yuanpeixiong/pangenome_heritability/test/kmers_directory \
    --threads 7



cp -r /storage/yangjianLab/yuanpeixiong/SV_project/sv_heritiability/test_data/alignment/temp_alignments/* /data/alignment/
