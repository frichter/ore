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
