#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Embryo data variant calling with GATK and samtools mpileup.

Felix Richter
felix.richter@icahn.mssm.edu
3/21/2019
Description: Run ORE through BSUB
"""


import os
import subprocess
import re
import glob

from ore_wrapper import OREwrapper


"""Set parameters."""
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
out_class = 'normal'  # rank normal extrema
out_prefix = home_dir + tissue + '_ore_100kb'  # _ore_100kb
# format order: sv_list out_class max_outs_list
# _outliers_5pct_max or just _outliers
outlier_output = (home_dir + tissue + '_outliers_5pct_max/' +
                  tissue + '_ore_SV{}_outliers_{}_lt{}.txt')
var_class_list = ['intronic', 'intergenic', 'exonic', 'UTR5',
                  'UTR3', 'splicing', 'upstream', 'ncRNA']
# allvars to ignore
# var_class = 'intergenic'  # 'allvars'  #
# exon_class_list = ['synonymous', 'nonsynonymous', '"frameshift|stopgain"']
# sv_list out_class max_outs_list var_class
# tss = '100 250 500 750 1000 2000 5000 1e4'
tss = '1e4 1e5'
enrich_f = (home_dir + tissue + '_enrich_heartenn_1e5/' + tissue +
            '_ens_ref_SV{}_{}_lt{}_{}_rmZ5pct_heartenn.txt')
# _rmZ5pct_CTCF_and_heart_TFs _enrich_map300
# atrial_enrich_map300_tss1e5

af_rare = '0.05 1e-2 1e-3 1e-4 1e-5 0.5'  #
af_min = '0 0 0 0 0 0.05'  #

rm_ids = '1-01013 1-01019 1-01094 1-02618 1-02702 1-04537 1-13670'

"""Annotations."""

heartenn_dir = ('/sc/orga/projects/chdiTrios/Felix/wgs/' +
                'pathogenicity_scores/deepheart_rare_var_scores/')
f_list = [i for i in glob.iglob(
    heartenn_dir + 'wgs_atrial_ids_max_score_filtered/*')]
f_list.sort()
# anno_str = ' '.join([i for i in glob.iglob(
#     heartenn_dir + 'wgs_atrial_ids_max_score_filtered/*')])
# annos = ('{0}wgs_atrial_ids_max_score.bed ' +
#          ' '.join(f_list[1:14]) + ' ').format(heartenn_dir)
# don't use all variants (for memory considerations)
annos = (' '.join(f_list[1:14]) + ' ').format(heartenn_dir)

"""Run ORE individually."""

# ore_obj = OREwrapper(home_dir, vcf, expr_f, out_class, out_prefix,
#                      outlier_output, var_class, enrich_f,
#                      # rm_ids: None for gtex, rm_ids for pcgc atrial
#                      rm_ids=rm_ids,
#                      tissue=tissue,
#                      annotations=annos,  # None annos
#                      tss=tss,
#                      exon_class=None,
#                      af_rare=af_rare, af_min=af_min)
# max_outs_i = ore_obj.max_outs_list[2]
# sv_i = '5'
# ore_cmd_w_args = ore_obj.run_ORE(sv_i, max_outs_i)
# ore_cmd_w_args = re.sub(' 2 2.5 3', ' 2', ore_cmd_w_args)
# # ore_cmd_w_args = re.sub('--processes 3', '--processes 1', ore_cmd_w_args)
# print(ore_cmd_w_args)  # + ' --n_perms 1000'  + ' --n_perms 1000'
# subprocess.call(ore_cmd_w_args, shell=True)

#

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
    ore_cmd_w_args = re.sub(' 2 2.5 3', ' 2', ore_cmd_w_args)
    subprocess.call(ore_cmd_w_args, shell=True)


#
