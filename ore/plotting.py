#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Plotting functions.

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-03-26
:Copyright: 2018, Felix Richter
:License: CC BY-SA
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt


def plot_outs_per_id(expr_outlier_df, outs_per_id_file):
    """Plot and save the number of outliers per ID as TXT and histogram.

    Args:
        expr_outlier_df (:obj:`DataFrame`): Pandas dataframe containing
            all the expression outlier data
        outs_per_id_file (:obj:`str`): prefix for histogram and TXT files

    """
    # aggregate number of outliers per ID
    outs_per_id = expr_outlier_df[['blinded_id', 'expr_outlier']].groupby(
        'blinded_id').sum()
    outs_per_id.to_csv(outs_per_id_file + '.txt', sep="\t")
    # plotting:
    median_outs_per_id = np.median(outs_per_id.expr_outlier)
    fig = plt.figure()
    plt.hist(outs_per_id.expr_outlier, color='black', bins=30)
    plt.axvline(x=median_outs_per_id, color='blue')
    plt.xlabel("Outliers per Sample")
    plt.ylabel("Number of Samples")
    plt.savefig(outs_per_id_file + '.png')
    plt.close(fig)


"""Figure out plotting enrichmnt confidence intervals

# from here: https://stats.stackexchange.com/a/2233
# https://stats.stackexchange.com/questions/1405/statistical-test-for-difference-between-two-odds-ratios

# loop over every line of dataframe
from numpy import exp, log, reciprocal
or_list = []
ci_list = []
enrich_df = pd.read_table()
enrich_df_line_x = enrich_df.iloc[0:1, :].tolist
or = enrich_df_line_x[11]
or_se = sqrt(sum(reciprocal(enrich_df_line_x.iloc[3:7])))
ci = (exp(log(or) - 1.96*or_se), exp(log(or) + 1.96*or_se))
ci_list.append(ci)

# plotting from here:
# http://hamelg.blogspot.com/2015/11/python-for-data-analysis-part-23-point.html
plt.figure(figsize=(9,9))

plt.errorbar(x=np.arange(0.1, 25, 1),
             y=sample_means,
             yerr=[(top-bot)/2 for top,bot in intervals],
             fmt='o')

plt.hlines(xmin=0, xmax=25,
           y=43.0023,
           linewidth=2.0,
           color="red")

"""
