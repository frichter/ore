#BSUB -W 6:00
#BSUB -q expressalloc
#BUSB -n 12
#BSUB -R "rusage[mem=5000]"
#BSUB -P acc_chdiTrios
#BSUB -J vent_outliers
#BSUB -m mothra
#BSUB -o vent_outliers.stdout
#BSUB -e vent_outliers.stderr


##################### PCGC Vent ##########################

VCF="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_01/wgs_vent_ids.norm.vcf.gz"
EXPR_F="/sc/orga/projects/chdiTrios/Felix/rna/pcgc/expression_data_rpkm_cutoff/ns_vent/residual_expr_5_SVs_hg19.bed.gz"
OUT_PREFIX="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_04/wgs_vent_expressAlloc"
ENRICH_PREFIX="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_04/enrichment_results/wgs_vent"

# vent_outliers


##################### PCGC Arterial/valve ##########################
# 
# VCF="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_01/wgs_arterial_valve_ids.norm.vcf.gz"
# EXPR_F="/sc/orga/projects/chdiTrios/Felix/rna/pcgc/expression_data_rpkm_cutoff/ns_art_valve_da/residual_expr_5_SVs_hg19.bed.gz"
# OUT_PREFIX="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_04/wgs_arterial_valve"
# ENRICH_PREFIX="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_04/enrichment_results/wgs_arterial_valve"

# art_valve_outliers



module load bedtools/2.27.0
module load samtools/1.3
module load bcftools/1.6
module load python/3.5.0
module load py_packages/3.5

cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore

# upstream and downstream (together) all annovar or subset
ENRICH_F=$ENRICH_PREFIX"/tssBi_SV5_norm_ncRNA_genes.txt"
python -m ore.ore --vcf $VCF \
    --bed $EXPR_F \
    --output $OUT_PREFIX \
    --outlier_output "outliers_norm_SV5.txt" \
    --enrich_file $ENRICH_F \
    --distribution "normal" \
    --threshold 2 3 4 \
    --max_outliers_per_id 1000 \
    --af_rare 0.05 1e-2 1e-3 1e-4 1e-5 \
    --tss_dist 1e4 \
    --annovar \
    --variant_class "ncRNA" \
    --humandb_dir "/sc/orga/projects/chdiTrios/whole_genome/humandb" \
    --processes 12

