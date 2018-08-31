#BSUB -W 6:00
#BSUB -q expressalloc
#BUSB -n 12
#BSUB -R "rusage[mem=10000]"
#BSUB -P acc_chdiTrios
#BSUB -J atrial_outs
#BSUB -m mothra
#BSUB -o atrial_outs.stdout
#BSUB -e atrial_outs.stderr


cd /sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_07

module load bedtools/2.27.0
module load samtools/1.3
module load bcftools/1.6
module load python/3.5.0
module load py_packages/3.5


PARENT_DIR="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_07"
VCF="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_01/wgs_atrial_ids.norm.vcf.gz"
EXPR_F="/sc/orga/projects/chdiTrios/Felix/rna/pcgc/expression_data_rpkm_cutoff/ns_atrial/residual_expr_5_SVs_hg19.bed.gz"
OUT_PREFIX="$PARENT_DIR/atrial_ore"
OUTLIER_OUT="$PARENT_DIR/atrial_ore_SV5_outliers.txt"
ENRICH_F="$PARENT_DIR/atrial_enrich/atrial_by_deepheart_anno_first100_10kb.txt"
# ore_per_anno_

# try these annotations:
# /sc/orga/projects/chdiTrios/Felix/wgs/bed_annotations
# ucsc_2017_03/factorbookMotif/CTCF.sorted.bed
# ucsc_2017_03/factorbookMotif/CTCF-ext.sorted.bed
# hg19_all/Centipedehg19.bed

TF_DIR="/sc/orga/projects/chdiTrios/Felix/wgs/bed_annotations/ucsc_2017_03"
ANNO_LIST="$TF_DIR/factorbookMotif/YY1.sorted.bed $TF_DIR/factorbookMotif/CTCF.sorted.bed $TF_DIR/TfbsClustered_split/CTCF.sorted.bed $TF_DIR/TfbsClustered_split/EP300.sorted.bed $TF_DIR/TfbsClustered_split/YY1.sorted.bed $TF_DIR/Conserved_TF_sites/NKX25.sorted.bed $TF_DIR/Conserved_TF_sites/P300.sorted.bed $TF_DIR/Conserved_TF_sites/YY1.sorted.bed"

cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore

# confirm on correct branch:
git status | head -n1
python -m ore.ore --version

# upstream and downstream (together) all annovar or subset
# python -m ore.ore --vcf $VCF \
time mprof run --include-children --multiprocess python -m ore.ore --vcf $VCF \
    --bed $EXPR_F \
    --output $OUT_PREFIX \
    --outlier_output $OUTLIER_OUT \
    --enrich_file $ENRICH_F \
    --distribution "normal" \
    --extrema \
    --threshold 2 \
    --max_outliers_per_id 1000 \
    --af_rare 0.05 1e-2 1e-3 1e-4 1e-5 \
    --intracohort_rare_ac 5 \
    --annotations $ANNO_LIST \
    --tss_dist 1e4 \
    --annovar \
    --humandb_dir "/sc/orga/projects/chdiTrios/whole_genome/humandb" \
    --processes 5


# --variant_class "UTR5" \
# --ensgene \
# --refgene \

mv atrial_ore_all_data.txt atrial_ore_all_data_extrema_deepheart_first100anno_10kb.txt
mv atrial_ore_rv_w_outliers.txt atrial_ore_rv_w_outliers_extrema_deepheart_first100anno_10kb.txt 

