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
module load python/3.6.2
module load py_packages/3.6

cd /sc/orga/projects/chdiTrios/Felix/dna_rna/rare_var_outliers/gtex_2018_08/

######################### GTEX LV #########################

EXPR_F="/hpc/users/richtf01/whole_genome/rare_variants_eqtl/gtex_control/gtex_final_expr_matrix/LV_gtex_2018_02_20/residual_expr_5_SVs_hg19.bed.gz"
VCF="/sc/orga/projects/chdiTrios/Felix/dna_rna/rare_var_outliers/gtex_june_2017/wgs_gtex.vcf.gz"
OUT_PREFIX="/sc/orga/projects/chdiTrios/Felix/dna_rna/rare_var_outliers/gtex_2018_08/lv_gtex"
ENRICH_PREFIX="/sc/orga/projects/chdiTrios/Felix/dna_rna/rare_var_outliers/gtex_2018_04/lv_gtex_enrich_utr5"

cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore

# upstream and downstream (together) all annovar or subset
ENRICH_F=$ENRICH_PREFIX"/tssBi_SV5_norm_ncRNA_refgene_ensgene.txt"
time python -m ore.ore --vcf $VCF \
    --bed $EXPR_F \
    --output $OUT_PREFIX \
    --outlier_output "outliers_norm_SV5.txt" \
    --enrich_file $ENRICH_F \
    --distribution "normal" \
    --threshold 2 \
    --max_outliers_per_id 500 \
    --af_rare 0.05 1e-2 1e-3 1e-4 1e-5 \
    --tss_dist 1e4 \
    --annovar \
    --variant_class "UTR5" \
    --ensgene \
    --refgene \
    --humandb_dir "/sc/orga/projects/chdiTrios/whole_genome/humandb" \
    --processes 6



### Custom virtual environment for ORE
module purge
module load python/3.5.0 py_packages/3.5
virtualenv venv_ore
source venv_ore/bin/activate
# confirm correct version
python --version
pip install statsmodels
deactivate

module purge
module load bedtools/2.27.0 samtools/1.3 bcftools/1.6
module load python/3.5.0 py_packages/3.5
source ~/venv_ore/bin/activate
## confirm 3.5
python --version

## did not work:
# module purge
# module load bedtools/2.27.0
# module load samtools/1.3
# module load bcftools/1.6
# module load python/3.6.2
# module load py_packages/3.6
# pip install memory_profiler

