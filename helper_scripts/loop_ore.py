"""Loop over ORE.

cd /sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_09

cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore

module purge
module load bedtools/2.27.0 samtools/1.3 bcftools/1.6
module load python/3.5.0 py_packages/3.5
source ~/venv_ore/bin/activate
## confirm 3.5
python --version
python

"""

import os
import subprocess
import re
import glob

from helper_scripts.ore_wrapper import OREwrapper


# change directories here
# cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore
os.chdir('/sc/orga/projects/chdiTrios/Felix/dna_rna/ore')

# atrial vent art_valve_da
tissue = 'atrial'
# wgs_pcgc_2018_09 for vent and art_valve_da, wgs_pcgc_2018_08 for atrial
if tissue in ['vent', 'art_valve_da']:
    home_dir = '/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_09/'
else:
    home_dir = '/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_08/'


vcf = (home_dir + '../wgs_pcgc_2018_01/wgs_' + tissue +
       '_ids.norm_smaller.vcf.gz')
expr_f = ('/sc/orga/projects/chdiTrios/Felix/rna/pcgc/expression_data_rpkm' +
          '_cutoff/ns_' + tissue + '/residual_expr_{}_SVs_hg19.bed.gz')
out_class = 'extrema'  # rank normal extrema
out_prefix = home_dir + tissue + '_ore_100kb'  # _ore_100kb
# format order: sv_list out_class max_outs_list
# _outliers_5pct_max or just _outliers
outlier_output = (home_dir + tissue + '_outliers_5pct_max/' +
                  tissue + '_ore_SV{}_outliers_{}_lt{}.txt')
var_class_list = ['intronic', 'intergenic', 'exonic', 'UTR5',
                  'UTR3', 'splicing', 'upstream', 'ncRNA']
# allvars to ignore
var_class = 'intergenic'  # 'allvars'  #
exon_class_list = ['synonymous', 'nonsynonymous', '"frameshift|stopgain"']
# sv_list out_class max_outs_list var_class
# tss = '100 250 500 750 1000 2000 5000 1e4'
tss = '1e4 1e5'
enrich_f = (home_dir + tissue + '_enrich_heartenn_1e5/' + tissue +
            '_ens_ref_SV{}_{}_lt{}_{}_rmZ5pct_heartenn_hiAF.txt')
# _rmZ5pct_CTCF_and_heart_TFs _enrich_map300
# atrial_enrich_map300_tss1e5

af_rare = '0.05 1e-2 1e-3'  # ' 1e-4 1e-5 0.5'  #
af_min = '0 0 0'  # ' 0 0 0.05'  #

rm_ids = '1-01013 1-01019 1-01094 1-02618 1-02702 1-04537 1-13670'
# This does not include the IDs excluded for arterial and valve analyses..
'1-00543 1-01399 1-06239 1-06552'  # vascular/valve
'1-00070 1-00713 1-01553 1-03547 1-04947 1-05433 1-06181 1-08783'  # vent

anno_dir = '/sc/orga/projects/chdiTrios/Felix/wgs/bed_annotations/'
annos = ('{0}ucsc_2017_03/factorbookMotif/CTCF.sorted.bed ' +
         '{0}hg19_all/any_*.bed ' +
         '{0}hg19_all/Centipedehg19.bed ' +
         '{0}hg19_all/DNaseMasterMajority.sorted.bed ').format(anno_dir)

# moving to atrial_ore_per_chrom subdirectory:
# mv var_* anno_high_prior_heart/
# mv *anno.bed anno_high_prior_heart/
# mv *anno_list.txt anno_high_prior_heart/
annos = ('{0}ucsc_2017_03/factorbookMotif/CTCF.sorted.bed ' +
         '{0}hg19_all/Centipedehg19.bed ' +
         '{0}hg19_all/DNaseMasterMajority.sorted.bed ' +
         '{0}hg19_all/E104*_coreMarks_* ' +
         # '{0}hg19_all/E08*_coreMarks_* ' +
         # '{0}hg19_all/E09*_coreMarks_* ' +
         '{0}hg19_all/*ata4* ' +
         '{0}hg19_all/heart%252c*ctss* ' +
         '{0}hg19_all/hg19.cage_peak_phase1and2combined_ann.bed ' +
         '{0}hg19_all/*kx2* ' +
         '{0}hg19_all/all_predictions.bed ' +
         '{0}hg19_all/robust_enhancers.sorted.bed ' +
         '{0}hg19_all/permissive_enhancers.sorted.bed ' +
         '{0}hg19_all/*bx5* ' +
         '{0}hg19_all/*.all_lncRNA.sorted.bed ').format(anno_dir)


heartenn_dir = ('/sc/orga/projects/chdiTrios/Felix/wgs/' +
                'pathogenicity_scores/deepheart_rare_var_scores/')
f_list = [i for i in glob.iglob(
    heartenn_dir + 'wgs_atrial_ids_max_score_filtered/*')]
f_list.sort()
# anno_str = ' '.join([i for i in glob.iglob(
#     heartenn_dir + 'wgs_atrial_ids_max_score_filtered/*')])
annos = ('{0}wgs_atrial_ids_max_score.bed ' +
         ' '.join(f_list[1:14]) + ' ').format(heartenn_dir)
annos = (' '.join(f_list[1:14]) + ' ').format(heartenn_dir)


"""
######### GTEx #########
"""

tissue = 'lv_gtex'
home_dir = ('/sc/orga/projects/chdiTrios/Felix/dna_rna/' +
            'rare_var_outliers/gtex_2018_08/')
vcf = ('/sc/orga/projects/chdiTrios/Felix/dna_rna/rare_var_outliers/' +
       'gtex_june_2017/wgs_gtex.vcf.gz')
expr_f = ('/hpc/users/richtf01/whole_genome/rare_variants_eqtl/gtex_control/' +
          'gtex_final_expr_matrix/LV_gtex_2018_02_20/' +
          'residual_expr_{}_SVs_hg19.bed.gz')
out_prefix = home_dir + 'lv_gtex'
# format order: sv_list out_class max_outs_list
outlier_output = (home_dir + 'lv_gtex_outs/' +
                  'lv_ore_SV{}_outliers_{}_lt{}.txt')
enrich_f = (home_dir + 'lv_gtex_enrich/' +
            'lv_gtex_ens_ref_SV{}_{}_lt{}_{}_rmZ5pct_wCommonVars.txt')
# wComVar2 wCommonVars

"""
########################################################################
# Running just a single example (e.g., rank-based)
########################################################################
"""

ore_obj = OREwrapper(home_dir, vcf, expr_f, out_class, out_prefix,
                     outlier_output, var_class, enrich_f,
                     # rm_ids: None for gtex, rm_ids for pcgc atrial
                     rm_ids=rm_ids,
                     tissue=tissue,
                     annotations=annos,  # None annos
                     tss=tss,
                     exon_class=None,
                     af_rare=af_rare, af_min=af_min)
max_outs_i = ore_obj.max_outs_list[2]
sv_i = '5'
ore_cmd_w_args = ore_obj.run_ORE(sv_i, max_outs_i)
ore_cmd_w_args = re.sub(' 2 2.5 3', ' 2', ore_cmd_w_args)
# ore_cmd_w_args = re.sub('--processes 3', '--processes 1', ore_cmd_w_args)
print(ore_cmd_w_args)  # + ' --n_perms 1000'  + ' --n_perms 1000'
subprocess.call(ore_cmd_w_args, shell=True)
# subprocess.call(ore_cmd_w_args + ' --n_perms 1', shell=True)

# change the name for rapid iterations
new_data_f = re.sub('wCommonVars', 'wCommonVars1', ore_obj.enrich_f_i)
mv_cmd = 'mv {} {}'.format(ore_obj.enrich_f_i, new_data_f)
subprocess.call(re.sub('.txt', '_variant.txt', mv_cmd), shell=True)
subprocess.call(re.sub('.txt', '_gene.txt', mv_cmd), shell=True)


"""
########################################################################
# LOOP OVER Variant classes
########################################################################
"""

for var_class_i in var_class_list:
    print(var_class_i)
    ore_obj = OREwrapper(home_dir, vcf, expr_f, out_class, out_prefix,
                         outlier_output, var_class_i, enrich_f,
                         rm_ids=rm_ids,
                         tissue=tissue,
                         annotations=annos,  # None annos
                         tss=tss,
                         exon_class=None,
                         af_rare=af_rare, af_min=af_min)
    max_outs_i = ore_obj.max_outs_list[2]
    sv_i = '5'
    ore_cmd_w_args = ore_obj.run_ORE(sv_i, max_outs_i)
    print(ore_cmd_w_args)
    subprocess.call(ore_cmd_w_args, shell=True)
    # if possible use the same all_data.txt file for all variants
    # after running, move the data to a permanent home so it is not overwritten
    # not necessary since you now put all data in the enrich directory
    # var_class_i = 'allvars'
    # new_data_f = (home_dir + tissue + '_data/' + tissue + '_ore_all_data_' +
    #               'SV{}_{}_lt{}_{}_rmZ5pct.txt').format(
    #               sv_i, out_class, max_outs_i, var_class_i)
    # mv_cmd = ore_obj.clean_files_after_run(new_data_f)
    # print(mv_cmd)
    # subprocess.call(mv_cmd, shell=True)

# loop over empty dataframes
anno_f_loc = out_prefix + '_per_chrom/tmp_bed_*.anno.bed'
for anno_f in glob.iglob(anno_f_loc):
    # only delete the zero-length overlap files
    if os.stat(anno_f).st_size != 0:
        continue
    print('deleting', anno_f)
    os.remove(anno_f)
    f_for_del = re.sub('anno.bed$', 'anno_list.txt', anno_f)
    print('deleting', f_for_del)
    os.remove(f_for_del)
    f_for_del = re.sub('anno.bed$', 'txt', anno_f)
    f_for_del = re.sub('tmp_bed_', 'var_', f_for_del)
    print('deleting', f_for_del)
    os.remove(f_for_del)


"""
########################################################################
# LOOP OVER EXON classes
########################################################################
"""

var_class_i = 'exonic'
ore_obj = OREwrapper(home_dir, vcf, expr_f, out_class, out_prefix,
                     outlier_output, var_class_i, enrich_f, rm_ids,
                     tissue)
max_outs_i = ore_obj.max_outs_list[2]
sv_i = '5'
# Only for arterial/valvar/da:
# ore_cmd_w_args = ore_obj.run_ORE(sv_i, max_outs_i)
# print(ore_cmd_w_args)
# subprocess.call(ore_cmd_w_args, shell=True)

new_data_f = (home_dir + tissue + '_data/' + tissue + '_ore_all_data_' +
              'SV{}_{}_lt{}_{}_rmZ5pct.txt').format(
              sv_i, out_class, max_outs_i, var_class_i)
cp_cmd = 'cp {} {}'.format(
    new_data_f, ore_obj.out_prefix + '_all_data_extrema.txt')
print(cp_cmd)
subprocess.call(cp_cmd, shell=True)

for exon_class_i in exon_class_list:
    print(exon_class_i)
    ore_obj = OREwrapper(home_dir, vcf, expr_f, out_class, out_prefix,
                         outlier_output, var_class_i, enrich_f, rm_ids,
                         tissue, exon_class_i)
    max_outs_i = ore_obj.max_outs_list[2]
    sv_i = '5'
    ore_cmd_w_args = ore_obj.run_ORE(sv_i, max_outs_i)
    print(ore_cmd_w_args)
    subprocess.call(ore_cmd_w_args, shell=True)


"""
########################################################################
# LOOP OVER MAXIMUM OUTLIERS.
########################################################################
"""
ore_obj = OREwrapper(home_dir, vcf, expr_f, out_class, out_prefix,
                     outlier_output, var_class, enrich_f, rm_ids, tissue)

sv_i = '5'
for max_outs_i in ore_obj.max_outs_list:
    ore_cmd_w_args = ore_obj.run_ORE(sv_i, max_outs_i)
    print(ore_cmd_w_args)
    subprocess.call(ore_cmd_w_args, shell=True)
    # after running, move the data to a permanent home so it is not overwritten
    new_data_f = (home_dir + 'atrial_data/atrial_ore_all_data_' +
                  'SV{}_{}_lt{}_{}_rmZ5pct.txt').format(
                  sv_i, out_class, max_outs_i, var_class)
    mv_cmd = ore_obj.clean_files_after_run(new_data_f)
    print(mv_cmd)
    subprocess.call(mv_cmd, shell=True)


"""
########################################################################
# LOOP OVER SURROGATE VARIBLES.
########################################################################
"""
ore_obj = OREwrapper(home_dir, vcf, expr_f, out_class, out_prefix,
                     outlier_output, var_class, enrich_f, rm_ids)

sv_list = [str(i) for i in range(4, 6)]  # (1, 11)

for sv_i in sv_list:
    max_outs_i = ore_obj.max_outs_list[2]
    ore_cmd_w_args = ore_obj.run_ORE(sv_i, max_outs_i)
    print(ore_cmd_w_args)
    subprocess.call(ore_cmd_w_args, shell=True)
    # after running, move the data to a permanent home so it is not overwritten
    new_data_f = (home_dir + 'atrial_data/atrial_ore_all_data_' +
                  'SV{}_{}_lt{}_{}_rmZ5pct.txt').format(
                  sv_i, out_class, max_outs_i, var_class)
    mv_cmd = ore_obj.clean_files_after_run(new_data_f)
    print(mv_cmd)
    subprocess.call(mv_cmd, shell=True)


# Other custom commands
sv_i = sv_list[5]
max_outs_i = ore_obj.max_outs_list[2]
ore_cmd_w_args = ore_obj.run_ORE(sv_i, max_outs_i)
print(ore_cmd_w_args)
new_data_f = (home_dir + 'atrial_data/atrial_ore_all_data_' +
              'SV{}_{}_lt{}_{}_rmZ5pct.txt').format(
              sv_i, out_class, max_outs_i, var_class)
mv_cmd = ore_obj.clean_files_after_run(new_data_f)

subprocess.call(
    'rm ' + re.sub('.txt', '_gene.txt', ore_obj.enrich_f_i), shell=True)
subprocess.call(
    'rm ' + re.sub('.txt', '_variant.txt', ore_obj.enrich_f_i), shell=True)
subprocess.call('rm ' + home_dir + ore_obj.full_data_f, shell=True)


"""
######### Alzheimer's Disease #########
cd /sc/orga/projects/chdiTrios/Felix/alzheimers/ore_2018_05
"""

tissue_dict = {'10': '19', '22': '20', '36': '17', '44': '17'}
tissue = '22'  # 10 22 36 44
sv_i = tissue_dict[tissue]
home_dir = ('/sc/orga/projects/chdiTrios/Felix/alzheimers/' +
            'ore_2018_05/')
vcf = home_dir + '../wgs/ad_wgs_cp.vcf.gz'
expr_f = (home_dir + '../expression/residuals_AMPAD_MSSM_GE_SV_{}' +
          '_tissue_' + tissue + '_with_disease_in_model_europeans_only_' +
          'wgs_ids_new_z.bed.gz')
# expr_f = (home_dir + '../expression/rnaseq_w_wgs_ids/residuals_AMPAD_' +
#           'MSSM_GE_SV_{}_tissue_' + tissue + '_with_disease_in_model_' +
#           'europeans_only_wgs_ids.bed.gz')
out_class = 'extrema'  # rank normal extrema
var_class_list = ['exonic',
                  'UTR3', 'splicing', 'upstream', 'ncRNA',
                  'UTR5',
                  'intronic', 'intergenic']
# allvars to ignore
var_class = 'UTR5'  # 'allvars'
# exon_class_list = ['synonymous', 'nonsynonymous', '"frameshift|stopgain"']
out_prefix = home_dir + 'ad_ore'
# format order: sv_list out_class max_outs_list
outlier_output = (home_dir + 'tissue_' + tissue + '_outs/' +
                  'SV{}_outliers_{}_lt{}.txt')
enrich_f = (home_dir + 'tissue_' + tissue + '_enrich/' +
            'ens_ref_SV{}_{}_lt{}_{}_rmZ5pct.txt')

rm_ids = None

# Single run
ore_obj = OREwrapper(home_dir, vcf, expr_f, out_class, out_prefix,
                     outlier_output, var_class, enrich_f, rm_ids,
                     tissue)
max_outs_i = '500'  # ore_obj.max_outs_list[2]
ore_cmd_w_args = ore_obj.run_ORE(sv_i, max_outs_i)
print(ore_cmd_w_args)  # + ' --n_perms 1000'  + ' --n_perms 1000'
subprocess.call(ore_cmd_w_args, shell=True)


"""
######### AD multitissue #########
cd /sc/orga/projects/chdiTrios/Felix/alzheimers/ore_2018_05
"""

home_dir = ('/sc/orga/projects/chdiTrios/Felix/alzheimers/' +
            'ore_2018_05/')
vcf = home_dir + '../wgs/ad_wgs_cp.vcf.gz'

# min2 min3 min4
tissue = 'min3'
sv_i = tissue
expr_f = (home_dir + '../expression/multitissue/' +
          'residuals_multitissue_' + tissue + 'tis_newZ_sorted.bed.gz')
out_class = 'extrema'  # rank normal extrema
var_class_list = ['exonic', 'UTR3', 'splicing', 'upstream', 'ncRNA',
                  'UTR5', 'intronic', 'intergenic']
# allvars to ignore
var_class = 'exonic'  # 'allvars'
# exon_class_list = ['synonymous', 'nonsynonymous', '"frameshift|stopgain"']
out_prefix = home_dir + 'ad_ore'
# format order: sv_list out_class max_outs_list
outlier_output = (home_dir + tissue + 'tissue_outs/' +
                  'SV{}_outliers_{}_lt{}.txt')
enrich_f = (home_dir + tissue + 'tissue_enrich/' +
            'ens_ref_SV{}_{}_lt{}_{}_rmZ5pct.txt')

rm_ids = None

# Single run
ore_obj = OREwrapper(home_dir, vcf, expr_f, out_class, out_prefix,
                     outlier_output, var_class, enrich_f, rm_ids,
                     tissue)
max_outs_i = '500'  # ore_obj.max_outs_list[2]
ore_cmd_w_args = ore_obj.run_ORE(sv_i, max_outs_i)
print(ore_cmd_w_args)  # + ' --n_perms 1000'  + ' --n_perms 1000'
subprocess.call(ore_cmd_w_args, shell=True)

"""
######### AD variant loop (single and multi tissue) #########
"""

for var_class_i in var_class_list:
    print(var_class_i)
    ore_obj = OREwrapper(home_dir, vcf, expr_f, out_class, out_prefix,
                         outlier_output, var_class_i, enrich_f, rm_ids,
                         tissue)
    max_outs_i = '500'  # ore_obj.max_outs_list[2]
    ore_cmd_w_args = ore_obj.run_ORE(sv_i, max_outs_i)
    if not os.path.exists(ore_obj.enrich_f_i[:-4] + '_gene.txt'):
        print(ore_cmd_w_args)
        subprocess.call(ore_cmd_w_args, shell=True)
    # new_data_f = (home_dir + 'tissue_' + tissue + '_data/' + tissue +
    #               '_ore_all_data_SV{}_{}_lt{}_{}_rmZ5pct.txt').format(
    #               sv_i, out_class, max_outs_i, var_class_i)
    # mv_cmd = ore_obj.clean_files_after_run(new_data_f)
    # print(mv_cmd)
    # subprocess.call(mv_cmd, shell=True)


"""
import subprocess
top_dir = ('/sc/orga/projects/chdiTrios/Felix/dna_rna/' +
           'wgs_pcgc_2018_08/atrial_ore_per_chrom')

chrom_list = [str(i) for i in range(1, 23)]
chrom_list = ['X', 'Y']
for chrom_i in chrom_list:
    cut_cmd = ('time cut -f-85,143-155,517- ' +
        '{top_dir}/all_tf_var_files/var_{chrom}.txt > ' +
        '{top_dir}/var_{chrom}.txt').format(
        top_dir=top_dir, chrom=chrom_i)
    print(cut_cmd)
    subprocess.call(cut_cmd, shell=True)

"""


"""
########################################################################
# Move permutations only if they do not exist
########################################################################

import subprocess
import os
import re
import glob

top_dir = ('/sc/orga/projects/chdiTrios/Felix/dna_rna/' +
           'wgs_pcgc_2018_09/atrial_ore_per_chrom')
perm_dir_list = ['perms2', 'perms3', 'perms4']  # 'perms',

for perm_dir in perm_dir_list:
    print(perm_dir)
    perm_iter = glob.iglob('{}/{}/*.txt'.format(top_dir, perm_dir))
    for perm_f in perm_iter:
        new_f = re.sub(perm_dir + '/', 'perms/', perm_f)
        if not os.path.exists(new_f):
            mv_cmd = 'mv {} {}'.format(perm_f, new_f)
            # print(mv_cmd)
            std_out = subprocess.call(mv_cmd, shell=True)
        else:
            print(new_f + ' already exists')

# cd /sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_08

## cross checking the factorbook results

anno_f = 'atrial_ore_all_data_extrema.txt'
tf_anno_f = 'atrial_ore_all_data_extrema_CTCF_heartTF.txt'
count = 0
with open(anno_f, 'r') as in_f, open(tf_anno_f, 'w') as out_f:
    header = next(in_f).strip().split('\t')
    for line in in_f:
        line_dict = dict(zip(header, line.strip().split('\t')))
        tf_list = [line_dict['factorbookMotif_CTCF'], line_dict['any_nkx25'],
                   line_dict['any_tbx5'], line_dict['any_gata4'],
                   line_dict['Centipedehg19'],
                   line_dict['DNaseMasterMajority']]
        if '1' in tf_list:
            _ = out_f.write(line)
        count += 1
        if count % 5000 == 0:
            print(count/20555356)


import pandas as pd
out_loc = ('atrial_outliers_5pct_max/' +
           'atrial_ore_SV5_outliers_extrema_lt200_forAnno.txt')
out_df = pd.read_table(out_loc)
out_genes = out_df[out_df.expr_outlier]['gene']

anno_df = pd.read_table(tf_anno_f)
anno_df.head()
anno_df.shape
min_intra_af = anno_df.intra_cohort_af.min()

af_cut_off = 1e-5
rv = (anno_df.popmax_af <= af_cut_off) & (
    anno_df.intra_cohort_af <= min_intra_af)
rv_df = anno_df[rv]

rv_df.head()
rv_df.shape

rv_outs_df = rv_df[rv_df.gene.isin(out_genes)]
rv_outs_df.head()
rv_outs_df.shape

fb_ctcf_rv = rv_outs_df.factorbookMotif_CTCF > 0
rv_outs_df[fb_ctcf_rv & rv_outs_df.expr_outlier]
fb_ctcf_rv.sum()
# cbind(c(2, 54), c(1717-2, (1717*81) - 54))

anno_rv = rv_outs_df.DNaseMasterMajority > 0
rv_outs_df[anno_rv & rv_outs_df.expr_outlier]
rv_outs_df[anno_rv][['gene', 'blinded_id']].drop_duplicates().shape

# gata4:
cbind(c(4, 110), c(1717-4, (1717*81) - 110))

# tbx5:
cbind(c(6, 150-6), c(1717-6, (1717*81) - 144))

# DNaseMasterMajority:
cbind(c(5, 140-5), c(1717-6, (1717*81) - 135))

# trouble shooting
import pandas as pd
f_loc = ('/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_08/' +
         'atrial_enrich_heartenn/all_data/atrial_ens_ref_SV5_extrema_' +
         'lt200_allvars_rmZ5pct_heartenn_all_data.txt')

df = pd.read_table(f_loc)
df.head()

sum(df.expr_outlier)
gt25neg1_vars = df.wgs_atrial_ids_max_score_gt25neg2 == 1
df_gt25neg1 = df[gt25neg1_vars]
rare_vars = df_gt25neg1.popmax_af <= 1e-5
df_gt25neg1[rare_vars].head()
df_gt25neg1.columns.values
df_gt25neg1[rare_vars].intra_cohort_af.min
sum(df_gt25neg1[rare_vars].intra_cohort_af <= 0.005618)
sum(rare_vars)  # should be 112+5

df.shape
df_not_scored = df[df.wgs_atrial_ids_max_score == 0]
df_not_scored.shape

"""

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
