#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Multi-tissue outlier calcluations

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-07-07
:Copyright: 2018, Felix Richter
:License: CC BY-SA

cd /sc/orga/projects/chdiTrios/Felix/alzheimers/expression
module load samtools/1.3
module load python/3.5.0
module load py_packages/3.5
python

"""

import glob
import re
# import functools
import pandas as pd
import subprocess as sp

expr_f_loc = "residuals_AMPAD*_only_wgs_ids_new_z.bed.gz"
expr_f_list = [i for i in glob.iglob(expr_f_loc)]
expr_f_list

expr_df_list = []
for expr_f in expr_f_list:
    # expr_f_list_i = expr_f_list[0]
    current_tissue = re.sub(".*_tissue_|_with_disease_in.*", "", expr_f)
    print(current_tissue)
    expr_df = pd.read_table(expr_f, sep='\t')
    print(expr_df.shape)
    expr_df.rename(columns={expr_df.columns[3]: "gene"}, inplace=True)
    expr_long = pd.melt(expr_df, id_vars=['gene', '#Chr', 'start', 'end'],
                        var_name='id', value_name='z_expr_' + current_tissue)
    # expr_long.set_index(['gene', '#Chr', 'start', 'end', 'id'], inplace=True)
    expr_df_list.append(expr_long)


# Trying sequential joins
for expr_df in expr_df_list:
    print(expr_df.head())
    expr_df.set_index(['gene', 'id'], inplace=True)
    # expr_df.set_index(['gene', '#Chr', 'start', 'end', 'id'], inplace=True)
    print(expr_df.head())
    expr_df.drop(['#Chr', 'start', 'end'], axis='columns', inplace=True)
    print(expr_df.head())
    expr_df.reset_index(inplace=True)

expr_joined = expr_df_list[0].join(expr_df_list[1], how='outer')
expr_joined.head()
expr_joined.shape
expr_joined = expr_joined.join(expr_df_list[2], how='outer')
expr_joined.head()
expr_joined.shape
expr_joined = expr_joined.join(expr_df_list[3], how='outer')
expr_joined.head()
expr_joined.shape

#
# expr_df_joined = functools.reduce(
#     lambda x, y: pd.merge(x, y, how='outer'), expr_df_list)
# for expr_df in expr_df_list:
#     print(expr_df.head())

# https://stackoverflow.com/a/33750531
# expr_joined['median'] = expr_joined.median(axis=1)
expr_joined = expr_joined.assign(median_z=expr_joined.median(axis=1))
# how many samples have all 4 tissues? 3, 2, 1?
expr_joined = expr_joined.assign(n_tissue=4 - expr_joined.isnull().sum(axis=1))
expr_long_loc = 'all_tissue_long_expr.txt'
expr_joined.to_csv(expr_long_loc, index=True, sep="\t", float_format='%g')

# Get gene info
expr_f = expr_f_list[0]
gene_info_df = pd.read_table(
    expr_f, sep='\t', usecols=[0, 1, 2, 3], index_col='gene')
gene_info_df.head()
expr_joined.reset_index(inplace=True)
expr_joined.set_index(['gene'], inplace=True)
expr_joined.head()
expr_joined.shape
expr_joined = gene_info_df.join(expr_joined, how='inner')
expr_joined.reset_index(inplace=True)

# expr_joined.rename(columns={'median': 'median_z'}, inplace=True)

min_tis = 4
expr_multi_tissues = expr_joined[expr_joined.n_tissue >= min_tis].drop(
    ['#Chr', 'start', 'end', 'z_expr_44', 'z_expr_22', 'z_expr_36',
     'z_expr_10', 'n_tissue'], axis='columns')
expr_multi_tissues.head()
expr_multi_tissues_wide = expr_multi_tissues.pivot_table(
    values='median_z', index='gene', columns='id')
expr_multi_tissues_wide.head()
expr_multi_tissues_wide.shape
expr_multi_tissues_wide = expr_multi_tissues_wide.join(
    gene_info_df, how='inner')
expr_multi_tissues_wide.reset_index(inplace=True)
cols = expr_multi_tissues_wide.columns.tolist()
cols
cols = cols[-3:] + cols[:-3]
expr_multi_tissues_wide = expr_multi_tissues_wide[cols]
multi_tissue_loc = 'residuals_multitissue_min{}tis_newZ.bed'.format(min_tis)
expr_multi_tissues_wide.to_csv(
    multi_tissue_loc, index=False, sep="\t", float_format='%g')

# sort, bgzip, and index
# expr_joined.set_index(['gene'], inplace=True)
min_tissue_list = [3, 4]
bed_f = 'residuals_multitissue_min{0}tis_newZ'
new_f_cmd = 'head -n1 {0}.bed > {0}_sorted.bed'
sort_cmd = ('sed 1d {0}.bed | sort -V -k1,1 -k2,2 >> {0}_sorted.bed')
zip_index_cmd = ('time bgzip {0}_sorted.bed && time tabix -p bed ' +
                 '{0}_sorted.bed.gz')
for min_tissue_i in min_tissue_list:
    print(min_tissue_i)
    bed_f = 'residuals_multitissue_min{0}tis_newZ'.format(min_tissue_i)
    print(new_f_cmd.format(bed_f))
    print(sort_cmd.format(bed_f))
    print(zip_index_cmd.format(bed_f))
    sp.call(new_f_cmd.format(bed_f), shell=True)
    sp.call(sort_cmd.format(bed_f), shell=True)
    sp.call(zip_index_cmd.format(bed_f), shell=True)


#
#
#
#
