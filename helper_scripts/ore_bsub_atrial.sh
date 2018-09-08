#BSUB -W 6:00
#BSUB -q alloc
#BUSB -n 12
#BSUB -R "rusage[mem=5000]"
#BSUB -P acc_chdiTrios
#BSUB -J atrial_outs
#BSUB -m mothra
#BSUB -o atrial_outs.stdout
#BSUB -e atrial_outs.stderr



### Creating a custom virtual environment for ORE
# module purge
# module load python/3.5.0 py_packages/3.5
# virtualenv venv_ore
# source venv_ore/bin/activate
# # confirm correct version
# python --version
# pip install statsmodels
# deactivate


cd /sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_09

module purge
module load bedtools/2.27.0 samtools/1.3 bcftools/1.6
module load python/3.5.0 py_packages/3.5
source ~/venv_ore/bin/activate
## confirm 3.5
python --version


# try these annotations:
TF_DIR="/sc/orga/projects/chdiTrios/Felix/wgs/bed_annotations/ucsc_2017_03"
# ANNO_LIST="$TF_DIR/factorbookMotif/YY1.sorted.bed $TF_DIR/factorbookMotif/CTCF.sorted.bed $TF_DIR/TfbsClustered_split/CTCF.sorted.bed $TF_DIR/TfbsClustered_split/EP300.sorted.bed $TF_DIR/TfbsClustered_split/YY1.sorted.bed $TF_DIR/Conserved_TF_sites/NKX25.sorted.bed $TF_DIR/Conserved_TF_sites/P300.sorted.bed $TF_DIR/Conserved_TF_sites/YY1.sorted.bed"
ANNO_LIST="$TF_DIR/factorbookMotif/*.sorted.bed"
 # $TF_DIR/TfbsClustered_split/*.sorted.bed $TF_DIR/Conserved_TF_sites/*.sorted.bed
###

## for annotations:
PARENT_DIR="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_08"

## for re-runs:
PARENT_DIR="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_09"
VCF="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_01/wgs_atrial_ids.norm_smaller.vcf.gz"
SV="10"
EXPR_F="/sc/orga/projects/chdiTrios/Felix/rna/pcgc/expression_data_rpkm_cutoff/ns_atrial/residual_expr_${SV}_SVs_hg19.bed.gz"
OUT_PREFIX="$PARENT_DIR/atrial_ore"
OUTLIER_OUT="$PARENT_DIR/atrial_ore_SV${SV}_outliers_norm_lt500.txt"
# atrial_ore_SV5_outliers_norm_lt500.txt
## removes 7 IDs (below) that are also removed for direct comparisons
# atrial_ore_SV5_outliers_extrema_customIDrm.txt
# atrial_ore_SV5_outliers_rank_customIDrm.txt
VAR_CLASS="UTR5"
ENRICH_F="$PARENT_DIR/atrial_enrich/atrial_enrich_norm_${VAR_CLASS}_SV${SV}_lt500.txt"
# ore_per_anno_
RM_IDS="1-01013 1-01019 1-01094 1-02618 1-02702 1-04537 1-13670"

##

cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore

# confirm on correct branch:
# git status | head -n1
# python -m ore.ore --help
# python -m ore.ore --version

# upstream and downstream (together) all annovar or subset
# time mprof run --include-children --multiprocess python -m ore.ore --vcf $VCF \
time python -m ore.ore --vcf $VCF \
    --bed $EXPR_F \
    --output $OUT_PREFIX \
    --outlier_output $OUTLIER_OUT \
    --enrich_file $ENRICH_F \
    --distribution "normal" \
    --threshold 2 \
    --max_outliers_per_id 500 \
    --af_rare 0.05 1e-2 1e-3 1e-4 1e-5 \
    --tss_dist 1e4 \
    --annovar \
    --variant_class "$VAR_CLASS" \
    --ensgene \
    --refgene \
    --humandb_dir "/sc/orga/projects/chdiTrios/whole_genome/humandb" \
    --processes 6




## profile: mprofile_20180906211602.dat mprofile_20180908115906.dat
real    1022m33.976s
user    7261m36.688s
sys     368m34.895s

# start: 2018-09-06 21:16:02
# 2018-09-07 14:18:35

## manually excluding IDs is faster
## for z-score, use --max_outliers_per_id 500
# --annotations $ANNO_LIST \
# --variant_class "UTR5" \
# --ensgene \
# --refgene \
# --extrema \
# time mprof run --include-children --multiprocess python -m ore.ore --vcf $VCF \
# --af_rare 0.05 1e-2 1e-3 1e-4 1e-5 \
# --annotations $ANNO_LIST \
# --threshold 2 \
# --intracohort_rare_ac 5 \
#     --threshold 0.025 0.01 \
# --max_outliers_per_id 500 \
# --exon_class nonsynonymous,synonymous,nonframeshift,frameshift,stopgain,stoploss
# --exon_class "nonsynonymous" \ "frameshift|stopgain"
## note that frameshift is first so that it's interpreted as ^frameshift
# --exclude_ids $RM_IDS \


echo $OUTLIER_OUT

cd /sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_09

## norm
mv atrial_ore_all_data.txt atrial_data/atrial_ore_all_data_lt500_SV${SV}_${VAR_CLASS}.txt
mv atrial_ore_rv_w_outliers.txt atrial_data/atrial_ore_rv_w_outliers_lt500_SV${SV}_${VAR_CLASS}.txt 

# most extreme
mv atrial_ore_all_data_extrema.txt atrial_ore_all_data_extrema_customIDrm.txt
mv atrial_ore_rv_w_outliers_extrema.txt atrial_ore_rv_w_outliers_extrema_customIDrm.txt 

#rank based
mv atrial_ore_all_data_rank.txt atrial_ore_all_data_rank_customIDrm.txt
mv atrial_ore_rv_w_outliers_rank.txt atrial_ore_rv_w_outliers_rank_customIDrm.txt 

## cleaning outliers
mv atrial_ore_small_vcf*outliers* atrial_outlier_data/


#
