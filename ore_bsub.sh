#BSUB -W 30:00
#BSUB -q alloc
#BUSB -n 12
#BSUB -R "rusage[mem=10000]"
#BSUB -P acc_chdiTrios
#BSUB -J gtex_lv
#BSUB -m mothra
#BSUB -o gtex_lv.stdout
#BSUB -e gtex_lv.stderr

module load bedtools/2.27.0
module load samtools/1.3
module load bcftools/1.6
module load python/3.5.0
module load py_packages/3.5

######################### AD Tissue 36 #########################

# cd /sc/orga/projects/chdiTrios/Felix/alzheimers

EXPR_F="/sc/orga/projects/chdiTrios/Felix/alzheimers/expression/residuals_AMPAD_MSSM_GE_SV_17_tissue_36_with_disease_in_model_europeans_only.bed.gz"
VCF_DIR="/sc/orga/projects/AMPADWGS/RawDataSinai/SCH_11923_06_16_2017/Project_SCH_11923_B02_GRM_WGS.2017-05-17/jgwd/joint_vcf/"
VCF=$VCF_DIR"SCH_11923_B02_GRM_WGS_2017-05-15_1.recalibrated_variants.vcf.gz"
OUT_PREFIX="/sc/orga/projects/chdiTrios/Felix/alzheimers/ore_2018_05/ad_ore"
ENRICH_F="/sc/orga/projects/chdiTrios/Felix/alzheimers/ore_2018_05/ad_ore_chr1_enrich_test.txt"

cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore

# upstream and downstream (together) all annovar or subset
ENRICH_F=$ENRICH_PREFIX"/tssBi_SV5_norm_ncRNA_refgene_ensgene.txt"
python -m ore.ore --vcf $VCF \
    --bed $EXPR_F \
    --output $OUT_PREFIX \
    --enrich_file $ENRICH_F \
    --distribution "normal" \
    --threshold 2 3 4 \
    --max_outliers_per_id 1000 \
    --af_rare 0.05 1e-2 1e-3 1e-4 1e-5 \
    --tss_dist 1e4 \
    --annovar \
    --variant_class "ncRNA" \
    --ensgene \
    --refgene \
    --humandb_dir "/sc/orga/projects/chdiTrios/whole_genome/humandb" \
    --processes 6



# cd /sc/orga/projects/chdiTrios/Felix/alzheimers

