#BSUB -W 2:00
#BSUB -q alloc
#BUSB -n 8
#BSUB -R "rusage[mem=30000]"
#BSUB -P acc_chdiTrios
#BSUB -J ad_ore_allvars
#BSUB -m mothra
#BSUB -o ad_ore_allvars.stdout
#BSUB -e ad_ore_allvars.stderr

cd /sc/orga/projects/chdiTrios/Felix/alzheimers/ore_2018_05

module load bedtools/2.27.0
module load samtools/1.3
module load bcftools/1.6
module load python/3.5.0
module load py_packages/3.5

######################### AD Tissue 36 #########################

# cd /sc/orga/projects/chdiTrios/Felix/alzheimers

PARENT_DIR="/sc/orga/projects/chdiTrios/Felix/alzheimers"
EXPR_F="$PARENT_DIR/expression/residuals_AMPAD_MSSM_GE_SV_17_tissue_36_with_disease_in_model_europeans_only_new_z_wgs_ids.bed.gz"
VCF="$PARENT_DIR/wgs/ad_wgs_cp.vcf.gz"
# ore_2018_06 ore_2018_05
OUT_PREFIX="$PARENT_DIR/ore_2018_05/ad_ore"
OUTLIER_OUT="$PARENT_DIR/ore_2018_05/most_extreme_outs_t36/ad_ore_outliers.txt"
ENRICH_F="$PARENT_DIR/ore_2018_05/most_extreme_t36_enrich/ad_ore_exonic_10kb_check_CIS.txt"

cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore

# upstream and downstream (together) all annovar or subset
# ENRICH_F=$ENRICH_PREFIX"/tssBi_SV5_norm_ncRNA_refgene_ensgene.txt"
python -m ore.ore --version

time mprof run --include-children --multiprocess python -m ore.ore --vcf $VCF \
    --bed $EXPR_F \
    --output $OUT_PREFIX \
    --outlier_output $OUTLIER_OUT \
    --enrich_file $ENRICH_F \
    --distribution "normal" \
    --extrema \
    --threshold 2 2.5 3 4 \
    --af_rare 5e-2 1e-2 1e-3 1e-4 1e-5 \
    --intracohort_rare_ac 5 \
    --tss_dist 1e3 2e3 5e3 1e4 \
    --annovar \
    --refgene \
    --ensgene \
    --variant_class "exonic" \
    --humandb_dir "/sc/orga/projects/chdiTrios/whole_genome/humandb" \
    --processes 5


#    --refgene \
#     --ensgene \
#     --variant_class "exonic" \
--exon_class "synonymous" \

mv ad_ore_all_data.txt ad_ore_all_data_exonic_ref_ens_10kb.txt
mv ad_ore_rv_w_outliers.txt ad_ore_rv_w_outliers_exonic_ref_ens_10kb.txt 

cp ad_ore_all_data_exonic_ref_ens_10kb.txt ad_ore_all_data.txt



# cd /sc/orga/projects/chdiTrios/Felix/alzheimers

