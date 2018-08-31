#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Recalculate Z-scores for every RNAseq sample

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
import subprocess

import pandas as pd

expr_f_list = [i for i in glob.iglob("residuals_AMPAD*_only_wgs_ids.bed.gz")]

expr_f_list_i = expr_f_list[3]
expr_df = pd.read_table(expr_f_list_i, sep='\t')
expr_df.shape
expr_df.head()
expr_df.rename(columns={expr_df.columns[3]: "gene"}, inplace=True)
expr_df.set_index(['gene'], inplace=True)

gene_df = expr_df.iloc[:, :3]
gene_df.head()
expr_df = expr_df.iloc[:, 3:]

expr_df.shape
expr_df_mean = expr_df.mean(axis=1)
expr_df_std = expr_df.std(axis=1)
expr_df = expr_df.sub(expr_df_mean, axis=0)
expr_df = expr_df.div(expr_df_std, axis=0)

expr_df.head()
expr_df.shape

expr_df = gene_df.join(expr_df, how='inner')
expr_df.head()
expr_df.shape

# reorder columns
expr_df.reset_index(inplace=True)
cols = expr_df.columns.tolist()
cols = cols[1:4] + cols[0:1] + cols[4:]
expr_df = expr_df.loc[:, cols]


expr_df_loc = re.sub(".bed.gz", "_new_z.bed", expr_f_list_i)
expr_df.to_csv(expr_df_loc, index=False, sep="\t", float_format='%g')

tbx_cmd = "time bgzip {} && time tabix -p bed {}.gz".format(
    expr_df_loc, expr_df_loc)
print(tbx_cmd)
subprocess.call(tbx_cmd, shell=True)
