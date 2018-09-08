#BSUB -W 30:00
#BSUB -q alloc
#BUSB -n 12
#BSUB -R "rusage[mem=10000]"
#BSUB -P acc_chdiTrios
#BSUB -J gtex_lv
#BSUB -m mothra
#BSUB -o gtex_lv.stdout
#BSUB -e gtex_lv.stderr

module purge
module load bedtools/2.27.0 samtools/1.3 bcftools/1.6
module load python/3.5.0 py_packages/3.5
source ~/venv_ore/bin/activate
## confirm 3.5
python --version


cd /sc/orga/projects/chdiTrios/Felix/dna_rna/rare_var_outliers/gtex_2018_08/

######################### GTEX LV #########################

EXPR_F="/hpc/users/richtf01/whole_genome/rare_variants_eqtl/gtex_control/gtex_final_expr_matrix/LV_gtex_2018_02_20/residual_expr_5_SVs_hg19.bed.gz"
VCF="/sc/orga/projects/chdiTrios/Felix/dna_rna/rare_var_outliers/gtex_june_2017/wgs_gtex.vcf.gz"
PARENT_DIR="/sc/orga/projects/chdiTrios/Felix/dna_rna/rare_var_outliers/gtex_2018_08"
OUT_PREFIX="$PARENT_DIR/lv_gtex"
ENRICH_F="$PARENT_DIR/lv_gtex_enrich_ref_ens_SV5_rank_noRM_enrich_utr5.txt"
OUTLIER_OUT="$PARENT_DIR/lv_gtex_ore_SV5_outliers_rank_lt500.txt"
# atrial_ore_SV5_outliers_norm_lt500.txt
# atrial_ore_SV5_outliers_rank_lt500.txt

cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore

# time mprof run --include-children --multiprocess python -m ore.ore --vcf $VCF \
time python -m ore.ore --vcf $VCF \
    --bed $EXPR_F \
    --output $OUT_PREFIX \
    --outlier_output $OUTLIER_OUT \
    --enrich_file $ENRICH_F \
    --distribution "rank" \
    --threshold 0.025 \
    --af_rare 0.05 1e-2 1e-3 1e-4 1e-5 \
    --tss_dist 1e4 \
    --annovar \
    --variant_class "UTR5" \
    --ensgene \
    --refgene \
    --humandb_dir "/sc/orga/projects/chdiTrios/whole_genome/humandb" \
    --processes 6


deactivate

#     --max_outliers_per_id 500 \

## profile: mprofile_20180906210843.dat, mprofile_20180907061724.dat, mprofile_20180907070745.dat, mprofile_20180907071553.dat
# location:
# /sc/orga/projects/chdiTrios/Felix/dna_rna/ore/profiles/gtex_lv_2018_09_07
# 
# real    445m36.571s
# user    3093m11.665s
# sys     245m11.485s
# 
# real    5m41.741s
# user    32m27.714s
# sys     2m0.015s
# 
# real    4m20.510s
# user    0m45.035s
# sys     0m31.207s
# 
# real    4m30.396s
# user    6m48.814s
# sys     0m44.187s
## view individually with mprof plot mprofile_20180907070745.dat

## start time: 2018-09-06 21:08:43
## Long formatting done: 2018-09-07 04:15:59

