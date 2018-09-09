#BSUB -W 6:00
#BSUB -q expressalloc
#BUSB -n 12
#BSUB -R "rusage[mem=10000]"
#BSUB -P acc_chdiTrios
#BSUB -J art_valve_da_outliers
#BSUB -m mothra
#BSUB -o art_valve_da_outliers.stdout
#BSUB -e art_valve_da_outliers.stderr


##################### PCGC Arterial/valve/da ##########################

# vent_outliers_3
cd /sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_09

module purge
module load bedtools/2.27.0 samtools/1.3 bcftools/1.6
module load python/3.5.0 py_packages/3.5
source ~/venv_ore/bin/activate
## confirm 3.5
python --version


PARENT_DIR="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_09"
VCF="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_01/wgs_arterial_valve_ids.norm_smaller.vcf.gz"
EXPR_F="/sc/orga/projects/chdiTrios/Felix/rna/pcgc/expression_data_rpkm_cutoff/ns_art_valve_da/residual_expr_5_SVs_hg19.bed.gz"
OUT_PREFIX="$PARENT_DIR/art_valve_da_ore_small_vcf"
OUTLIER_OUT="$PARENT_DIR/art_valve_da_ore_small_vcf_SV5_outliers_norm_lt500.txt"
VAR_CLASS="exonic"
ENRICH_F="$PARENT_DIR/art_valve_da_enrich/art_valve_da_enrich_norm_${VAR_CLASS}_SV5_lt500.txt"

cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore

# confirm on correct branch:
# git checkout heart_ore_v27plus
# git status | head -n1
# python -m ore.ore --version

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
    --processes 3


deactivate


cd /sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_09

VAR_CLASS="allvars"

mv art_valve_da_ore_small_vcf_all_data.txt art_valve_da_data/art_valve_da_ore_small_vcf_all_data_lt500_SV5_${VAR_CLASS}.txt
mv art_valve_da_ore_small_vcf_rv_w_outliers.txt art_valve_da_data/art_valve_da_ore_rv_w_outliers_small_vcf_lt500_SV5_${VAR_CLASS}.txt 


# mprofile_20180908120510.dat mprofile_20180908162852.dat mprofile_20180908183531.dat
# real    248m6.220s
# user    1803m38.323s
# sys     184m1.028s
# 
# real    0m26.493s
# user    2m13.892s
# sys     0m16.565s
# 
# real    26m54.535s
# user    42m3.656s
# sys     4m4.846s
# 248+26+((6+26+54)/60)
# 275

# variant abstraction step time:
2018-09-08 12:05:13
2018-09-08 16:05:33
