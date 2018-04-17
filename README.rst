OutAR
-----

To use, confirm the following are installed:
    - module load bedtools/2.27.0
    - module load samtools/1.3
    - module load bcftools/1.6
    - module load python/3.5.0


Example run
    outar
    --vcf test.vcf.gz \
    --bed test.bed.gz \
    --enrich_file enrichment.txt \
    --distribution "normal" \
    --threshold 2 \
    --max_outliers_per_id 500 \
    --af_rare 0.05 \
    --tss_dist 5000

