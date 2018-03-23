#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Outlier calling.

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-02-08
:Copyright: 2018, Felix Richter
:License: CC BY-SA
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from genes import RNASeqError


def applyParallel(dfGrouped, func):
    """Parallelize the pandas apply function.

    SOURCE
        https://stackoverflow.com/questions/26187759/parallelize-apply-after-pandas-groupby

    TODO
        handle ValueError: No objects to concatenate
    """
    from multiprocessing import Pool, cpu_count
    print("Using all {} cores".format(cpu_count()))
    with Pool(cpu_count()) as p:
        ret_list = p.map(func, [group for name, group in dfGrouped])
    return pd.concat(ret_list)


class Outliers(object):
    """Methods and attributes of outliers."""

    def __init__(self, gene_pheno_loc):
        """Initialize outlier dataframe.

        Convert gene expression data frame from wide to long.
        """
        gene_expr_df = pd.read_table(gene_pheno_loc, low_memory=False)
        gene_expr_df = gene_expr_df.iloc[:, 3:]
        self.expr_long_df = pd.melt(gene_expr_df, id_vars='gene',
                                    value_vars=gene_expr_df.columns[1:],
                                    var_name='blinded_id', value_name='z_expr')

    def get_outliers(self, ids_to_keep, distribution="normal",
                     most_extreme=False, expr_cut_off=2):
        """Calculate RNAseq outliers.

        Options for dist=["normal", "rank", "custom"].

        Loads RNAseq in same format as for FastQTL, then identifies outliers.

        TODO:
            Consider making a static method
            Clean/incorporate parallelization method and/or ask stackoverflow
                on options for vectorizing the `find_expr_outlier` function
            Confirm expr_cut_off is an acceptable value for normal/rank/custom
            Check that WGS IDs are consistent with outlier IDs
        """
        self.expr_long_df = self.expr_long_df[
            self.expr_long_df.blinded_id.isin(ids_to_keep)]
        if distribution == "normal":
            if not self.test_normality():
                print("Not normal, so calculating z-scores...")
            print("Calculating z-score outliers....")
            self.identify_outliers_from_normal(expr_cut_off)
        elif distribution == "rank":
            print("Calculating ranks...")
            self.expr_long_df = applyParallel(self.expr_long_df.groupby(
                'gene'), self.calculate_ranks)
            print("Ranks calculated, identifying outliers")
            self.expr_long_df = self.identify_outliers_from_ranks(
                self.expr_long_df, expr_cut_off=0.05)
            print("Outliers identified, saving...")
        elif distribution == "custom":
            not_0_1 = ~self.expr_long_df.z_expr.isin([0, 1])
            if any(not_0_1):
                print(self.expr_long_df[not_0_1])
                raise RNASeqError("The values above were not 0 or 1")
            self.expr_long_df["expr_outlier"] = self.expr_long_df.z_expr == 1
        else:
            raise RNASeqError("'{}' is not a valid outlier distribution".
                              format(distribution))
        if most_extreme:
            print("Identifying most extreme outlier per gene...")
            expr_outlier_df = applyParallel(self.expr_long_df.groupby(
                'gene'), self.find_most_extreme_expr_outlier)
            expr_outlier_df['expr_outlier_neg'] = (
                expr_outlier_df.expr_outlier_neg &
                expr_outlier_df.expr_outlier)
        else:
            expr_outlier_df = self.expr_long_df
        expr_outlier_df.reset_index(inplace=True)
        return expr_outlier_df

    def identify_outliers_from_normal(self, expr_cut_off):
        """Identify outliers more extreme than a z-score threshold."""
        self.expr_long_df["z_abs"] = abs(self.expr_long_df.z_expr)
        self.expr_long_df["expr_outlier"] = (
            self.expr_long_df.z_abs > expr_cut_off)
        self.expr_long_df["expr_outlier_neg"] = (
            (self.expr_long_df.z_expr < 0) &
            self.expr_long_df.expr_outlier)

    @staticmethod
    def identify_outliers_from_ranks(expr_long_df, expr_cut_off):
        """Identify outliers based on those more extreme than percentile.

        Args
            `expr_cut_off`; percentile cut-off for outliers

        TODO
            confirm minimum expr_cut_off has at least 1 sample past
            the threshold
            Calculate outliers as those at least 3 IQRs from median

        """
        min_expr_cut_off = min(set(expr_long_df.expr_rank))
        if expr_cut_off <= min_expr_cut_off or expr_cut_off >= 0.5:
            raise RNASeqError("The percentile cut-off specified ({}) is " +
                              "not between 0.5 and the minimum cut-off " +
                              "for this sample size, {}".format(
                                expr_cut_off, min_expr_cut_off))
        hi_expr_cut_off = 1 - expr_cut_off
        expr_long_df["expr_outlier_neg"] = (
            expr_long_df.expr_rank <= expr_cut_off)
        expr_long_df["expr_outlier"] = (
            (expr_long_df.expr_rank >= hi_expr_cut_off) |
            expr_long_df.expr_outlier_neg)
        return expr_long_df

    @staticmethod
    def calculate_ranks(gene_group):
        """Calculate ranks for each gene.

        Args
            `gene_group`
        Returns
            `gene_group`: with expr_rank which is the percentile

        """
        gene_group["expr_rank"] = gene_group["z_expr"].rank(method='average',
                                                            pct=True)
        return gene_group

    def test_normality(self):
        """Check if each gene has normal distribution.

        Options include QQ-plots, shapiro-wilk test and others.

        """
        return True

    @staticmethod
    def find_most_extreme_expr_outlier(gene_group):
        """Label outliers in a group.

        Create a column for absolute value of expression z-score, determine
            which `blinded_id` has the maximum expression z-score to
            create `expr_outlier_status` column, then determine if any of these
            outliers have z-scores < 0 (i.e., are low or negative outliers).

        Args:
            `gene_group`: long-format expression dataframe for a gene

        TODO:
            Only applies to normal distribution, either generalize to ranks
            Or specify another way of identify max/min ranked per gene
        """
        gene_group['expr_outlier'] = (
            (gene_group.z_abs == max(gene_group.z_abs)) &
            gene_group.expr_outlier)
        # gene_group['expr_outlier_neg'] = (
        #     gene_group.expr_outlier_neg & gene_group.expr_outlier)
        return gene_group

    @staticmethod
    def get_ids_w_low_out_ct(expr_outlier_df, outlier_max):
        """Identify blinded_ids with a ton of outliers.

        Args:
            `expr_outlier_df`: long-format expression dataframe
                labeling each gene-ID as an outlier
            `outlier_max`: maximum number of outliers per ID

        """
        outs_per_id = expr_outlier_df[['blinded_id', 'expr_outlier']].groupby(
            'blinded_id').sum()
        ids_w_hi_out_ct = list(
            outs_per_id[outs_per_id.expr_outlier >= outlier_max].index)
        print("The following IDs have >= {} outliers each: {}".format(
            outlier_max, ", ".join(ids_w_hi_out_ct)))
        ids_to_keep = list(
            outs_per_id[outs_per_id.expr_outlier < outlier_max].index)
        return ids_to_keep

    @staticmethod
    def plot_out_per_id(expr_outlier_df, outs_per_id_file):
        """Plot and save the number of outliers per ID as TXT and histogram."""
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


#
