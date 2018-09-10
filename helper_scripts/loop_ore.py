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

from helper_scripts.ore_wrapper import OREwrapper


# atrial vent art_valve_da
tissue = 'art_valve_da'
home_dir = '/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_09/'
vcf = (home_dir + '../wgs_pcgc_2018_01/wgs_' + tissue +
       '_ids.norm_smaller.vcf.gz')
expr_f = ('/sc/orga/projects/chdiTrios/Felix/rna/pcgc/expression_data_rpkm' +
          '_cutoff/ns_' + tissue + '/residual_expr_{}_SVs_hg19.bed.gz')
out_class = 'extrema'  # rank normal extrema
out_prefix = home_dir + tissue + '_ore'
# sv_list out_class max_outs_list
outlier_output = (home_dir + tissue + '_outliers_5pct_max/' +
                  tissue + '_ore_SV{}_outliers_{}_lt{}.txt')
var_class_list = ['allvars', 'intronic', 'intergenic', 'exonic', 'UTR5',
                  'UTR3', 'splicing', 'upstream', 'ncRNA']
var_class = 'UTR5'
# sv_list out_class max_outs_list var_class
enrich_f = (home_dir + tissue + '_enrich_map300/' + tissue +
            '_ens_ref_SV{}_{}_lt{}_{}_rmZ5pct.txt')

rm_ids = '1-01013 1-01019 1-01094 1-02618 1-02702 1-04537 1-13670'


# change directories here
# cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore
os.chdir('/sc/orga/projects/chdiTrios/Felix/dna_rna/ore')

#


"""LOOP OVER Variant classes"""
var_class_i = var_class_list[0]
ore_obj = OREwrapper(home_dir, vcf, expr_f, out_class, out_prefix,
                     outlier_output, var_class_i, enrich_f, rm_ids,
                     tissue)
max_outs_i = ore_obj.max_outs_list[2]
sv_i = '5'
ore_cmd_w_args = ore_obj.run_ORE(sv_i, max_outs_i)
print(ore_cmd_w_args)
subprocess.call(ore_cmd_w_args, shell=True)

for var_class_i in var_class_list:
    ore_obj = OREwrapper(home_dir, vcf, expr_f, out_class, out_prefix,
                         outlier_output, var_class_i, enrich_f, rm_ids)
    max_outs_i = ore_obj.max_outs_list[2]
    sv_i = '5'
    ore_cmd_w_args = ore_obj.run_ORE(sv_i, max_outs_i)
    print(ore_cmd_w_args)
    subprocess.call(ore_cmd_w_args, shell=True)
    # should use the same all_data.txt file for all variants

# after running, move the data to a permanent home so it is not overwritten
new_data_f = (home_dir + tissue + '_data/' + tissue + '_ore_all_data_' +
              'SV{}_{}_lt{}_{}_rmZ5pct.txt').format(
              sv_i, out_class, max_outs_i, 'all_vars')
mv_cmd = ore_obj.clean_files_after_run(new_data_f)
print(mv_cmd)
subprocess.call(mv_cmd, shell=True)

"""LOOP OVER MAXIMUM OUTLIERS."""


"""LOOP OVER SURROGATE VARIBLES."""
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
