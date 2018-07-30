#BSUB -W 6:00
#BSUB -q expressalloc
#BUSB -n 12
#BSUB -R "rusage[mem=10000]"
#BSUB -P acc_chdiTrios
#BSUB -J vent_outliers_3
#BSUB -m mothra
#BSUB -o vent_outliers_3.stdout
#BSUB -e vent_outliers_3.stderr

######################### GTEX LV #########################

# EXPR_F="/hpc/users/richtf01/whole_genome/rare_variants_eqtl/gtex_control/gtex_final_expr_matrix/LV_gtex_2018_02_20/residual_expr_5_SVs_hg19.bed.gz"
# VCF="/sc/orga/projects/chdiTrios/Felix/dna_rna/rare_var_outliers/gtex_june_2017/wgs_gtex.vcf.gz"
# OUT_PREFIX="/sc/orga/projects/chdiTrios/Felix/dna_rna/rare_var_outliers/gtex_2018_04/wgs_gtex_lv_2"
# ENRICH_PREFIX="/sc/orga/projects/chdiTrios/Felix/dna_rna/rare_var_outliers/gtex_2018_04/enrichment_results/lv/"

# job name: gtex_lv_2


##################### PCGC Vent ##########################

# vent_outliers_3
cd /sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_07

module load bedtools/2.27.0
module load samtools/1.3
module load bcftools/1.6
module load python/3.5.0
module load py_packages/3.5


PARENT_DIR="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_07"
VCF="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_01/wgs_vent_ids.norm.vcf.gz"
EXPR_F="/sc/orga/projects/chdiTrios/Felix/rna/pcgc/expression_data_rpkm_cutoff/ns_vent/residual_expr_5_SVs_hg19.bed.gz"
OUT_PREFIX="$PARENT_DIR/vent_ore"
OUTLIER_OUT="$PARENT_DIR/vent_ore_SV5_utliers_most_extreme.txt"
ENRICH_F="$PARENT_DIR/vent_enrich_most_extreme/vent_ref_ens_splicing_10kb.txt"
# vent_ore_per_anno_10kb.txt

cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore

# confirm on correct branch:
# git checkout heart_ore_v27plus
git status | head -n1
python -m ore.ore --version

python -m ore.ore --vcf $VCF \
    --bed $EXPR_F \
    --output $OUT_PREFIX \
    --outlier_output $OUTLIER_OUT \
    --enrich_file $ENRICH_F \
    --distribution "normal" \
    --extrema \
    --threshold 2 2.5 3 4 \
    --max_outliers_per_id 1000 \
    --af_rare 0.05 1e-2 1e-3 1e-4 1e-5 \
    --intracohort_rare_ac 5 \
    --tss_dist 1e4 \
    --annovar \
    --humandb_dir "/sc/orga/projects/chdiTrios/whole_genome/humandb" \
    --processes 5



# --ensgene \
# --refgene \
# --variant_class "UTR5" \

mv vent_ore_all_data.txt vent_ore_all_data_utr5_ref_ens_10kb.txt
mv vent_ore_rv_w_outliers.txt vent_ore_rv_w_outliers_utr5_ref_ens_10kb.txt 


