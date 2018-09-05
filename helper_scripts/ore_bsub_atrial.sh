#BSUB -W 6:00
#BSUB -q alloc
#BUSB -n 12
#BSUB -R "rusage[mem=5000]"
#BSUB -P acc_chdiTrios
#BSUB -J atrial_outs
#BSUB -m mothra
#BSUB -o atrial_outs.stdout
#BSUB -e atrial_outs.stderr


cd /sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_08

module purge
module load bedtools/2.27.0
module load samtools/1.3
module load bcftools/1.6
# module load python/3.5.0 # no statsmodels
# module load py_packages/3.5
module load python/3.6.2 # no mprof
# statsmodels/0.8.0 which throws FutureWarning, so use statsmodels/0.9.0
module load py_packages/3.6


PARENT_DIR="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_08"
VCF="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_01/wgs_atrial_ids.norm.vcf.gz"
EXPR_F="/sc/orga/projects/chdiTrios/Felix/rna/pcgc/expression_data_rpkm_cutoff/ns_atrial/residual_expr_5_SVs_hg19.bed.gz"
OUT_PREFIX="$PARENT_DIR/atrial_ore"
OUTLIER_OUT="$PARENT_DIR/atrial_ore_SV5_outliers_extrema_customIDrm.txt"
# atrial_ore_SV5_outliers_norm_lt500.txt
## removes 7 IDs (below) that are also removed for direct comparisons
# atrial_ore_SV5_outliers_extrema_customIDrm.txt
# atrial_ore_SV5_outliers_rank_customIDrm.txt
ENRICH_F="$PARENT_DIR/atrial_enrich_ens_ref_norm_SV5.txt"
# ore_per_anno_
RM_IDS="1-01013 1-01019 1-01094 1-02618 1-02702 1-04537 1-13670"

# try these annotations:
# /sc/orga/projects/chdiTrios/Felix/wgs/bed_annotations
# ucsc_2017_03/factorbookMotif/CTCF.sorted.bed
# ucsc_2017_03/factorbookMotif/CTCF-ext.sorted.bed
# hg19_all/Centipedehg19.bed

TF_DIR="/sc/orga/projects/chdiTrios/Felix/wgs/bed_annotations/ucsc_2017_03"
ANNO_LIST="$TF_DIR/factorbookMotif/YY1.sorted.bed $TF_DIR/factorbookMotif/CTCF.sorted.bed $TF_DIR/TfbsClustered_split/CTCF.sorted.bed $TF_DIR/TfbsClustered_split/EP300.sorted.bed $TF_DIR/TfbsClustered_split/YY1.sorted.bed $TF_DIR/Conserved_TF_sites/NKX25.sorted.bed $TF_DIR/Conserved_TF_sites/P300.sorted.bed $TF_DIR/Conserved_TF_sites/YY1.sorted.bed"


##

cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore

# confirm on correct branch:
# git status | head -n1
# python -m ore.ore --help
# python -m ore.ore --version

# upstream and downstream (together) all annovar or subset
time python -m ore.ore --vcf $VCF \
    --bed $EXPR_F \
    --output $OUT_PREFIX \
    --outlier_output $OUTLIER_OUT \
    --enrich_file $ENRICH_F \
    --distribution "normal" \
    --threshold 2 \
    --exclude_ids $RM_IDS \
    --af_rare 0.05 1e-2 1e-3 1e-4 1e-5 \
    --tss_dist 1e4 \
    --annovar \
    --variant_class "upstream" \
    --ensgene \
    --refgene \
    --humandb_dir "/sc/orga/projects/chdiTrios/whole_genome/humandb" \
    --processes 3


## manually excluding IDs is faster
## for z-score, use --max_outliers_per_id 500
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



cd /sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_08

echo $OUTLIER_OUT

## norm
mv atrial_ore_all_data.txt atrial_ore_all_data_utr5_lt100.txt
mv atrial_ore_rv_w_outliers.txt atrial_ore_rv_w_outliers_utr5_lt100.txt 

# most extreme
mv atrial_ore_all_data_extrema.txt atrial_ore_all_data_extrema_utr5_customIDrm.txt
mv atrial_ore_rv_w_outliers_extrema.txt atrial_ore_rv_w_outliers_extrema_utr5_customIDrm.txt 

#rank based
mv atrial_ore_all_data_rank.txt atrial_ore_all_data_rank_utr5_customIDrm.txt
mv atrial_ore_rv_w_outliers_rank.txt atrial_ore_rv_w_outliers_rank_utr5_customIDrm.txt 



#
