#BSUB -W 60:00
#BSUB -q alloc
#BUSB -n 12
#BSUB -R "rusage[mem=10000]"
#BSUB -P acc_chdiTrios
#BSUB -J ad_ore
#BSUB -m mothra
#BSUB -o ad_ore.stdout
#BSUB -e ad_ore.stderr

module load bedtools/2.27.0
module load samtools/1.3
module load bcftools/1.6
module load python/3.5.0
module load py_packages/3.5

######################### AD Tissue 36 #########################

# cd /sc/orga/projects/chdiTrios/Felix/alzheimers

PARENT_DIR="/sc/orga/projects/chdiTrios/Felix/alzheimers"
EXPR_F="$PARENT_DIR/expression/residuals_AMPAD_MSSM_GE_SV_17_tissue_36_with_disease_in_model_europeans_only.bed.gz"
VCF="$PARENT_DIR/wgs/ad_wgs_cp.vcf.gz"
# ore_2018_06 ore_2018_05
OUT_PREFIX="$PARENT_DIR/ore_2018_06_job/ad_ore"
ENRICH_F="$PARENT_DIR/ore_2018_06_job/ad_ore_enrich_test.txt"

cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore

# upstream and downstream (together) all annovar or subset
# ENRICH_F=$ENRICH_PREFIX"/tssBi_SV5_norm_ncRNA_refgene_ensgene.txt"
python -m ore.ore --version

time mprof run --include-children --multiprocess python -m ore.ore --vcf $VCF \
    --bed $EXPR_F \
    --output $OUT_PREFIX \
    --enrich_file $ENRICH_F \
    --distribution "normal" \
    --extrema \
    --threshold 2 \
    --max_outliers_per_id 1000 \
    --af_rare 5e-2 1e-2 1e-3 \
    --tss_dist 5e4 \
    --annovar \
    --variant_class "UTR5" \
    --ensgene \
    --refgene \
    --humandb_dir "/sc/orga/projects/chdiTrios/whole_genome/humandb" \
    --processes 8



# cd /sc/orga/projects/chdiTrios/Felix/alzheimers

