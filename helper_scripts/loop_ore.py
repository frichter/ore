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


home_dir = '/sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_09/'
vcf = home_dir + '../wgs_pcgc_2018_01/wgs_atrial_ids.norm_smaller.vcf.gz'
sv_list = [str(i) for i in range(1, 11)]
expr_f = ('/sc/orga/projects/chdiTrios/Felix/rna/pcgc/expression_data_rpkm' +
          '_cutoff/ns_atrial/residual_expr_{}_SVs_hg19.bed.gz')
out_class = 'extrema'  # rank normal extrema
out_prefix = home_dir + 'atrial_ore'
# sv_list out_class max_outs_list
outlier_output = (home_dir + 'atrial_outliers_5pct_max/' +
                  'atrial_ore_SV{}_outliers_{}_lt{}.txt')
var_class = 'UTR5'
# sv_list out_class max_outs_list var_class
enrich_f = (home_dir + 'atrial_enrich_map300/' +
            'atrial_ens_ref_SV{}_{}_lt{}_{}_rmZ5pct.txt')

rm_ids = '1-01013 1-01019 1-01094 1-02618 1-02702 1-04537 1-13670'


# change directories here
# cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore
os.chdir('/sc/orga/projects/chdiTrios/Felix/dna_rna/ore')

#
ore_obj = OREwrapper(home_dir, vcf, expr_f, out_class, out_prefix,
                     outlier_output, var_class, enrich_f, rm_ids)

sv_i = sv_list[7]
max_outs_i = ore_obj.max_outs_list[2]
ore_cmd_w_args = ore_obj.run_ORE(sv_i, max_outs_i)
print(ore_cmd_w_args)
subprocess.call(ore_cmd_w_args, shell=True)

new_data_f = (home_dir + 'atrial_data/atrial_ore_all_data_' +
              'SV{}_{}_lt{}_{}_rmZ5pct.txt').format(
              sv_i, out_class, max_outs_i, var_class)
mv_cmd = ore_obj.clean_files_after_run(new_data_f)
print(mv_cmd)
subprocess.call(ore_cmd_w_args, shell=True)

# for norm:
if out_class is 'normal':
    max_outs_list = ['100', '250', '500', '750', '1000', '2500', '5000']
elif out_class is 'extrema':
    max_outs_list = ['100', '150', '200', '400', '600', '800', '1000']
else:
    max_outs_list = ['100', '150', '200', '400', '600', '800', '1000']


if out_class is 'extrema':
    extrema_arg = '--extrema '
    dist_arg = 'normal'
    max_outs_i = max_outs_list[2]
    # max_outs_arg = '--max_outliers_per_id ' + max_outs_i + ' '
    # rm_id_arg = ''
    """ Alternatively preset which IDs to remove:"""
    max_outs_i = 'custom'
    max_outs_arg = ''
    rm_id_arg = '--exclude_ids ' + rm_ids + ' '
elif out_class is 'normal':
    extrema_arg = ''
    dist_arg = out_class
    max_outs_i = max_outs_list[2]
    max_outs_arg = '--max_outliers_per_id ' + max_outs_i + ' '
    rm_id_arg = ''
    """ Alternatively preset which IDs to remove:"""
    # max_outs_i = 'custom'
    # max_outs_arg = ''
    # rm_id_arg = '--exclude_ids ' + rm_ids + ' '
elif out_class is 'rank':
    extrema_arg = ''
    dist_arg = out_class
    max_outs_i = 'custom'
    max_outs_arg = ''
    rm_id_arg = '--exclude_ids ' + rm_ids + ' '

expr_f_i = expr_f.format(sv_i)
outlier_output_i = outlier_output.format(sv_i, out_class, max_outs_i)
enrich_f_i = enrich_f.format(sv_i, out_class, max_outs_i, var_class)
ore_cmd = ('time python -m ore.ore --vcf {vcf} --bed {expr} ' +
           '--output {out_pref} --outlier_output {outlier_pref} ' +
           '--enrich_file {enrich} --distribution {dist} --threshold 2 ' +
           '{extrema_arg}{max_outs_arg}{rm_id_arg}' +
           '--af_rare 0.05 1e-2 1e-3 1e-4 1e-5 --tss_dist 5e3 1e4 ' +
           '--annovar --variant_class {var} --refgene --ensgene ' +
           '--humandb_dir "/sc/orga/projects/chdiTrios/whole_genome/' +
           'humandb" --processes 3')
ore_cmd_w_args = ore_cmd.format(
    vcf=vcf, expr=expr_f_i, out_pref=out_prefix,
    outlier_pref=outlier_output_i, enrich=enrich_f_i,
    dist=dist_arg, extrema_arg=extrema_arg,
    max_outs_arg=max_outs_arg, rm_id_arg=rm_id_arg,
    var=var_class)

print(ore_cmd_w_args)
subprocess.call(ore_cmd_w_args, shell=True)

"""Move intermediate files to a subdirectory so they can be remade."""

new_data_f = (home_dir + 'atrial_data/atrial_ore_all_data_' +
              'SV{}_{}_lt{}_{}_rmZ5pct.txt').format(
              sv_i, out_class, max_outs_i, var_class)
if out_class is 'extrema':
    full_data_f = 'atrial_ore_all_data_extrema.txt'
elif out_class is 'normal':
    full_data_f = 'atrial_ore_all_data.txt'
elif out_class is 'rank':
    full_data_f = 'atrial_ore_all_data_rank.txt'

mv_cmd = 'mv {} {}'.format(
    home_dir + full_data_f,
    new_data_f)
if not os.path.exists(new_data_f):
    print(mv_cmd)
    subprocess.call(mv_cmd, shell=True)


# Other custom commands
subprocess.call('rm ' + re.sub('.txt', '_gene.txt', enrich_f_i), shell=True)
subprocess.call('rm ' + re.sub('.txt', '_variant.txt', enrich_f_i), shell=True)
subprocess.call('rm ' + home_dir + full_data_f, shell=True)

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
