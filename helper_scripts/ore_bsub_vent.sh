#BSUB -W 6:00
#BSUB -q expressalloc
#BUSB -n 12
#BSUB -R "rusage[mem=10000]"
#BSUB -P acc_chdiTrios
#BSUB -J vent_outliers_3
#BSUB -m mothra
#BSUB -o vent_outliers_3.stdout
#BSUB -e vent_outliers_3.stderr


##################### PCGC Vent ##########################

# vent_outliers_3
cd /sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_09

module purge
module load bedtools/2.27.0 samtools/1.3 bcftools/1.6
module load python/3.5.0 py_packages/3.5
source ~/venv_ore/bin/activate
## confirm 3.5
python --version


PARENT_DIR="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_09"
VCF="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_01/wgs_vent_ids.norm_smaller.vcf.gz"
EXPR_F="/sc/orga/projects/chdiTrios/Felix/rna/pcgc/expression_data_rpkm_cutoff/ns_vent/residual_expr_5_SVs_hg19.bed.gz"
OUT_PREFIX="$PARENT_DIR/vent_ore_small_vcf"
OUTLIER_OUT="$PARENT_DIR/vent_ore_small_vcf_SV5_outliers_norm_lt500.txt"
ENRICH_F="$PARENT_DIR/vent_enrich_norm_utr5_SV5_lt500.txt"
# vent_ore_per_anno_10kb.txt

cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore

# confirm on correct branch:
# git checkout heart_ore_v27plus
# git status | head -n1
python -m ore.ore --version

# python -m ore.ore --vcf $VCF \
time mprof run --include-children --multiprocess python -m ore.ore --vcf $VCF \
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
    --variant_class "UTR5" \
    --ensgene \
    --refgene \
    --humandb_dir "/sc/orga/projects/chdiTrios/whole_genome/humandb" \
    --processes 12



# mprofile_20180907075150.dat mprofile_20180908103508.dat mprofile_20180908110016.dat
real    301m16.596s
user    1897m23.308s
sys     193m35.795s

real    22m31.307s
user    152m21.255s
sys     6m8.458s

real    4m34.601s
user    22m41.917s
sys     1m22.994s

# start: 2018-09-07 07:51:50
# VCF long done: 2018-09-07 12:53:02

# --ensgene \
# --refgene \
# --variant_class "UTR5" \

mv vent_ore_all_data.txt vent_ore_all_data_utr5_ref_ens_10kb.txt
mv vent_ore_rv_w_outliers.txt vent_ore_rv_w_outliers_utr5_ref_ens_10kb.txt 


