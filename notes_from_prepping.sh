


############################################################
# data preparation
############################################################

cut -f4-5 gene_TSS_0b_hg19.txt  > gene_strand_hg19.txt
cut -f4-5 gene_TSS_0b_hg38.txt  > gene_strand_hg38.txt
gzip gene_strand_hg*


http://annovar.openbioinformatics.org/en/latest/user-guide/download/
# registered and downloaded ANNOVAR data

perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
# perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad_genome humandb/
# last command taking a while

wget http://www.openbioinformatics.org/annovar/download/hg19_gnomad_genome.txt.idx.gz
wget http://www.openbioinformatics.org/annovar/download/hg19_gnomad_genome.txt.gz

time perl annovar/table_annovar.pl tests/test_data/test_per_chrom/tmp_annovar_input_chr1.txt
    -buildver hg19 \
    -out min_annovar_data_needed \
    --otherinfo -remove -protocol \
    refGene,ensGene,gnomad_genome \
    -operation g,g,f -nastring NA

## test with annovar on Minerva

############################################################
# preparing test files
############################################################

# also see tool_runs_2.sh and notes_for_testing.py
TESTDATADIR="/Users/felixrichter/Dropbox/PhD/ore/ore/tests/"
cd $TESTDATADIR
VCF_IN="test_data_all_atrial_samples/test.vcf.gz"
VCF_OUT="test_data/test.vcf.gz"
SAMPLES="/Users/felixrichter/Dropbox/PhD/rna_seq_other/expression_data_rpkm_cutoff/ns_atrial/t350_ids_only.txt"

# obtained from bcftools_commands.sh
time bcftools view --samples-file $SAMPLES --force-samples -Oz -o $VCF_OUT $VCF_IN
time tabix -p vcf $VCF_OUT

cp test_data_all_atrial_samples/test.bed.gz* test_data/

############################################################
# test runs, coverage, profiling, etc
############################################################

# subset_loc = "/sc/orga/projects/chdiTrios/Felix/dna_rna/outar_test.bed"

# TESTDATADIR="/sc/orga/projects/chdiTrios/Felix/dna_rna/outar_test/"
TESTDATADIR="/Users/felixrichter/Dropbox/PhD/ore/ore/tests/test_data"
VCF="$TESTDATADIR/test.vcf.gz"
EXPR_F="$TESTDATADIR/test.bed.gz"

cd /Users/felixrichter/Dropbox/PhD/ore
pip3 install -e .
cd $TESTDATADIR

# pip3 install ore

ore --version
ore --help

# or python3 -m ore.ore 

ore --vcf $VCF \
    --bed $EXPR_F \
    --enrich_file "$TESTDATADIR/test_enrichment.txt" \
    --distribution "normal" \
    --threshold 2 2.5 \
    --max_outliers_per_id 500 \
    --af_rare 0.25 0.1 0.05 \
    --tss_dist 2e4

 # rm -rf test_*


############################################################
# Preparing a readme/docs, documentation
############################################################

# https://stackoverflow.com/questions/40859607/add-usage-help-of-command-line-tool-to-readme-rst
ore --help > USAGE.rst
sed -e '/Usage/r USAGE.rst' -e '$G' README.template.rst > README2.rst
# confirm new changes are okay and edit README2 if not
diff README.rst README2.rst
# replace old w new:
mv README2.rst README.rst

# readthedocs
# http://docs.readthedocs.io/en/latest/getting_started.html#in-rst
cd docs
sphinx-quickstart # accept defaults for most

# convert to html (from within docs)
make html
# https://readthedocs.org/dashboard/ore/version/latest/

############################################################
# uploading to Pypi
############################################################

## final commands run: 
## need to update version EVERYTIME in version.py
# confirm readme rst is parsed correctly:
python3 setup.py check --restructuredtext
# confirm previous builds are archived
ls dist/ore*
mv dist/ore-0.1.* dist/archive/
# create new builds
python3 setup.py sdist
python3 setup.py bdist_wheel
twine upload dist/ore*


## previous commands run:
python3 setup.py sdist

# wheel is faster for installing that installing from the source
pip3 install wheel
# do not use --universal flag because only runs on python3 (ie not python2)
python3 setup.py bdist_wheel
# displays ttwo warnings:
# warning: no previously-included files found matching '.gitignore'
# warning: no previously-included files found matching '.coverage'

# deleted the following from setup.cfg
[bdist_wheel]
universal=1

# temporarily
mv /Users/felixrichter/.pypirc $HOME/.pypirc_temp 

# check the readme renders appropriately
# https://github.com/pypa/readme_renderer
pip3 install readme_renderer
python3 setup.py check -r -s

# practicing with https://test.pypi.org/

python3 setup.py sdist bdist_wheel

twine upload --repository-url https://test.pypi.org/legacy/ dist/*

twine upload --repository-url https://upload.pypi.org/legacy/ dist/*
twine upload --repository-url https://pypi.org/project/ore/ dist/*

twine register dist/*
twine upload dist/*



cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore
bsub < ore_bsub.sh

###############################################################
# Back on Minerva with real data. How to install python
# package locally? Just use as script
###############################################################

cd /sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_04/

module load bedtools/2.27.0
module load samtools/1.3
module load bcftools/1.6
module load python/3.5.0
module load py_packages/3.5

cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore

# bsub < ore_bsub.sh
# bsub < ore_bsub_atrial.sh 
# python -m ore.ore --help

# https://stackoverflow.com/questions/11536764/how-to-fix-attempted-relative-import-in-non-package-even-with-init-py

# upstream and downstream (together) all annovar or subset
# _ensgene_refgene
ENRICH_F=$ENRICH_PREFIX"/tssBi_SV5_norm_All_by_anno_test_new.txt"
time python -m ore.ore --vcf $VCF \
    --bed $EXPR_F \
    --output $OUT_PREFIX \
    --outlier_output "outliers_norm_SV5.txt" \
    --enrich_file $ENRICH_F \
    --distribution "normal" \
    --threshold 2 \
    --max_outliers_per_id 1000 \
    --af_rare 5e-2 1e-2 1e-3 \
    --tss_dist 1e4 \
    --annovar \
    --variant_class "UTR5" \
    --ensgene \
    --refgene \
    --humandb_dir "/sc/orga/projects/chdiTrios/whole_genome/humandb" \
    --processes 2


# --af_rare 0.05 1e-2 1e-3 1e-4 1e-5 \
# --variant_class "UTR5" \
# --ensgene \
# --refgene \
## append this command for profiling:
# time mprof run --include-children --multiprocess 
# e.g.:
# time mprof run --include-children --multiprocess python -m ore.ore --version

##################### PCGC Atrial ##########################

VCF="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_01/wgs_atrial_ids.norm.vcf.gz"
EXPR_F="/sc/orga/projects/chdiTrios/Felix/rna/pcgc/expression_data_rpkm_cutoff/ns_atrial/residual_expr_5_SVs_hg19.bed.gz"
OUT_PREFIX="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_04/wgs_atrial"
ENRICH_PREFIX="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_04/enrichment_results/wgs_atrial"


##################### PCGC Vent ##########################

VCF="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_01/wgs_vent_ids.norm.vcf.gz"
EXPR_F="/sc/orga/projects/chdiTrios/Felix/rna/pcgc/expression_data_rpkm_cutoff/ns_vent/residual_expr_5_SVs_hg19.bed.gz"
OUT_PREFIX="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_04/wgs_vent"
ENRICH_PREFIX="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_04/enrichment_results/wgs_vent"


##################### PCGC Arterial/valve ##########################

VCF="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_01/wgs_arterial_valve_ids.norm.vcf.gz"
EXPR_F="/sc/orga/projects/chdiTrios/Felix/rna/pcgc/expression_data_rpkm_cutoff/ns_art_valve_da/residual_expr_5_SVs_hg19.bed.gz"
OUT_PREFIX="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_04/wgs_arterial_valve"
ENRICH_PREFIX="/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_04/enrichment_results/wgs_arterial_valve"



####################
##### AD prep
####################


cd /sc/orga/projects/chdiTrios/Felix/alzheimers

module load bedtools/2.27.0
module load samtools/1.3
module load bcftools/1.6
module load python/3.5.0
module load py_packages/3.5

# concat from bcftools (https://samtools.github.io/bcftools/bcftools.html)
# into a single directory

ls /sc/orga/projects/AMPADWGS/RawDataSinai/SCH_11923_06_16_2017/Project_SCH_11923_B02_GRM_WGS.2017-05-17/jgwd/joint_vcf/*variants.vcf.gz
# copy and paste into text edit and manually sort into correct order
# then paste into vcf_files_by_chrom.txt in vim
# confirm sample headers are all exactly the same:
zgrep -m1 '#CHROM'  /sc/orga/projects/AMPADWGS/RawDataSinai/SCH_11923_06_16_2017/Project_SCH_11923_B02_GRM_WGS.2017-05-17/jgwd/joint_vcf/*variants.vcf.gz | sed 's/.*vcf.gz//g' > vcf_per_chrom_headers.txt
cat vcf_per_chrom_headers.txt | uniq | wc -l
# should be 1

cd /sc/orga/projects/chdiTrios/Felix/alzheimers/wgs
time bcftools concat -f vcf_files_by_chrom.txt --output ad_wgs.vcf.gz -Oz --threads 3
# user 586m50.518s
## get a warning, [W::hts_idx_load2] The index file is older than the data file: /sc/orga/projects/AMPADWGS/RawDataSinai/SCH_11923_06_16_2017/Project_SCH_11923_B02_GRM_WGS.2017-05-17/jgwd/joint_vcf/SCH_11923_B02_GRM_WGS_2017-05-15_others.recalibrated_variants.vcf.gz.tbi
cp --verbose /sc/orga/projects/AMPADWGS/RawDataSinai/SCH_11923_06_16_2017/Project_SCH_11923_B02_GRM_WGS.2017-05-17/jgwd/joint_vcf/*variants.vcf.gz .
## so copy and paste files to current directory, recreate tabix
for F in *_variants.vcf.gz; do time tabix -f -p vcf ${F}; done
# for loop from https://www.biostars.org/p/259000/
# ls *variants.vcf.gz
time bcftools concat -f vcf_files_by_chrom_cp.txt -Ou  --threads 5 | time bcftools view -f PASS -i "F_MISSING <= 0.3 && QUAL >= 30" -o ad_wgs_cp.vcf.gz -Oz --threads 5
# user 773m9.340s

cd /sc/orga/projects/chdiTrios/Felix/alzheimers/expression
for F in *.bed.gz; do time tabix -f -p bed ${F}; done

time tabix -p vcf ad_wgs_cp.vcf.gz

####################
##### AD RUN
####################


cd /sc/orga/projects/chdiTrios/Felix/alzheimers

module load bedtools/2.27.0
module load samtools/1.3
module load bcftools/1.6
module load python/3.5.0
module load py_packages/3.5

# memory_profiler package not available in py_packages/3.6
# module load python/3.6.2
# module load py_packages/3.6


cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore

python -m ore.ore --help
python -m ore.ore --version

PARENT_DIR="/sc/orga/projects/chdiTrios/Felix/alzheimers"
EXPR_F="$PARENT_DIR/expression/residuals_AMPAD_MSSM_GE_SV_17_tissue_36_with_disease_in_model_europeans_only.bed.gz"
VCF="$PARENT_DIR/wgs/ad_wgs_cp.vcf.gz"
# ore_2018_06 ore_2018_05
OUT_PREFIX="$PARENT_DIR/ore_2018_05/ad_ore"
ENRICH_F="$PARENT_DIR/ore_2018_05/ad_ore_enrich_test.txt"

time mprof run --include-children --multiprocess python -m ore.ore --vcf $VCF \
    --bed $EXPR_F \
    --output $OUT_PREFIX \
    --enrich_file $ENRICH_F \
    --distribution "normal" \
    --extrema \
    --threshold 2 \
    --max_outliers_per_id 1000 \
    --af_rare 5e-2 1e-2 1e-3 \
    --intracohort_rare_ac 5 \
    --tss_dist 5e4 \
    --annovar \
    --variant_class "UTR5" \
    --ensgene \
    --refgene \
    --humandb_dir "/sc/orga/projects/chdiTrios/whole_genome/humandb" \
    --processes 6

    # --outlier_output "outliers_norm_SV5.txt" \


# fixing up results from a previous run
rename "chr" "" tmp_*
time sed -i 's/chr//g' tmp_*

## run on node with internet
screen -R -D ore
module purge
module load python/3.5.0 py_packages/3.5
virtualenv venv_ore
source venv_ore/bin/activate
pip install --upgrade pip
pip install ore
deactivate


# Download ANNOVAR databases by (a) registering here:
http://www.openbioinformatics.org/annovar/annovar_download_form.php

# Will get a link like this with a specific USER_KEY
wget http://www.openbioinformatics.org/annovar/download/USER_KEY/annovar.latest.tar.gz
tar -xzvf annovar.latest.tar.gz

# Set your data directory. REPLACE WITH YOUR OWN human_db location
DB_DIR="./annovar/human_db/"

# This provides the custom ANNOVAR scripts
annotate_variation.pl
retrieve_seq_from_fasta.pl

DB_DIR="/sc/orga/projects/chdiTrios/whole_genome/humandb/"
cd annovar
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad_genome humandb/

## Run these 3 commands
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene $DB_DIR/
perl annotate_variation.pl --buildver hg19 --downdb seq $DB_DIR/hg19_seq
perl retrieve_seq_from_fasta.pl $DB_DIR/hg19_knownGene.txt -seqdir $DB_DIR/hg19_seq -format knownGene -outfile $DB_DIR/hg19_knownGeneMrna.fa


# Switch to node where you want to actually execute code
ssh interactive_node
screen -R -D ore
module purge
module load bedtools/2.27.0 samtools/1.3 bcftools/1.6
module load python/3.5.0 py_packages/3.5 
source venv_ore/bin/activate

# Run ORE!
VCF="wgs/ad_wgs_cp.vcf.gz"
OUTPUT_PREFIX="brain_region_1"
DB_DIR="annovar/human_db/"

ore --vcf $VCF \
    --bed $BED \
    --output $OUTPUT_PREFIX \
    --distribution "custom" \
    --af_rare 5e-2 1e-2 1e-3 1e-4 1e-5 \
    # --intracohort_rare_allele_count 5 \
    --tss_dist 5e3 \
    --annovar \
    --humandb_dir $DB_DIR \
    --n_processes 6


time tar -zcvf wgs_by_chrom.tar.gz wgs_by_chrom/
# VCF_DIR="/sc/orga/projects/AMPADWGS/RawDataSinai/SCH_11923_06_16_2017/Project_SCH_11923_B02_GRM_WGS.2017-05-17/jgwd/joint_vcf/"
# VCF=$VCF_DIR"SCH_11923_B02_GRM_WGS_2017-05-15_1.recalibrated_variants.vcf.gz"

