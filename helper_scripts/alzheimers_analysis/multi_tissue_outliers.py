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
import functools
import pandas as pd

expr_f_loc = "residuals_AMPAD*_only_wgs_ids_new_z.bed.gz"
expr_f_list = [i for i in glob.iglob(expr_f_loc)]

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


expr_df_joined = functools.reduce(
    lambda x, y: pd.merge(x, y, how='outer'), expr_df_list)
for expr_df in expr_df_list:
    print(expr_df.head())

# https://stackoverflow.com/a/33750531
expr_df_joined['median'] = expr_df_joined.median(axis=1)
expr_df_joined = expr_df_joined.assign(median=expr_df_joined.median(axis=1))

# how many samples have all 4 tissues? 3, 2, 1?

#
#
#
#
