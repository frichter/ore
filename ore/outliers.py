#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Outlier calling.

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-02-08
:Copyright: 2018, Felix Richter
:License: CC BY-SA
"""

import os
import re

import pandas as pd

from .genes import RNASeqError
from .plotting import plot_outs_per_id
from .utils import applyParallel


class Outliers(object):
    """Methods and attributes of outliers."""

    def __init__(self, pheno_loc, output_prefix, outlier_postfix,
                 extrema, distribution, threshold, logger):
        """Initialize outlier dataframe.

        Args:
            pheno_loc (:obj:`str`): gene expression (phenotype) location
            output_prefix (:obj:`str`): file prefix for outputs
            outlier_postfix (:obj:`str`): file ending for outlier files
            extrema (:obj:`boolean`): T/F for using most extreme outlier
            distribution (:obj:`str`): type of outlier distribution considered
            threshold (:obj:`list`): list of outlier cut-off thresholds
            logger (:obj:`logging object`): Current logger

        Attributes:
            expr_long_df (:obj:`DataFrame`): RNAseq expression in long format
            expr_outs_loc (:obj:`str`): full outlier file location
            extrema (:obj:`boolean`): T/F for using most extreme outlier
            distribution (:obj:`str`): type of outlier distribution considered
            threshold (:obj:`list`): list of outlier cut-off thresholds
            least_extr_threshold (:obj:`float`): least extreme threshold

        Raises:
            :obj:`RNASeqError`: if `distribution` is not valid, i.e., not in
                ["normal", "rank", "custom"]

        TODO:
            low_memory is deprecated so instead specify any ambiguous dtypes
            https://stackoverflow.com/questions/24251219/pandas-read-csv-low-memory-and-dtype-options

        """
        gene_expr_df = pd.read_table(pheno_loc, low_memory=False)
        logger.debug(gene_expr_df.head())
        logger.debug(gene_expr_df.columns)
        logger.debug(gene_expr_df.columns.values[2])
        logger.debug(gene_expr_df.shape)
        gene_expr_df = gene_expr_df.iloc[:, 3:]
        # Convert gene expression data frame from wide to long:
        self.expr_long_df = pd.melt(
            gene_expr_df,
            id_vars=gene_expr_df.columns.values[2],  # 'gene',
            value_vars=gene_expr_df.columns[1:].tolist(),
            var_name='blinded_id',
            value_name='z_expr')
        # set the output file location
        self.expr_outs_loc = (output_prefix + "_outliers.txt")
        if outlier_postfix:
            self.expr_outs_loc = (output_prefix + "_" + outlier_postfix)
        # set states on which specific outlier definitions are being used
        self.extrema = extrema
        self.distribution = distribution
        self.threshold = threshold
        if isinstance(self.threshold, float):
            self.least_extr_threshold = threshold
        elif self.distribution == "normal":
            self.least_extr_threshold = min(self.threshold)
        elif self.distribution == "rank":
            self.least_extr_threshold = max(self.threshold)
        elif self.distribution == "custom":
            self.least_extr_threshold = 1
        else:
            raise RNASeqError("'{}' is not a valid outlier distribution".
                              format(distribution))

    def prepare_outliers(self, outlier_max, vcf_id_list):
        """Obtain gene expression outliers.

        Args:
            outlier_max (:obj:`int`): maximum number of outliers per ID
            vcf_id_list (:obj:`list`): list of Blinded IDs with WGS

        Check if expression outlier file already exists. If it does not,
            use `get_outliers` to obtain outliers using only IDs
            with WGS. Once these outliers are obtained, plot a histogram
            of the number of outliers per ID. If there is a maximum
            number of outliers per ID specified, then remove IDs that
            cross this threshold and call outliers again. Save the
            outlier dataframe to `expr_outs_loc`

        """
        if os.path.exists(self.expr_outs_loc):
            return "Already made outlier file " + self.expr_outs_loc
        # only work with IDs in WGS that are also in the RNAseq
        lines_w_consistent_ids = self.expr_long_df.blinded_id.isin(vcf_id_list)
        if lines_w_consistent_ids.shape[0] == 0:
            raise RNASeqError("No overlapping IDs between RNAseq and VCF")
        self.expr_long_df = self.expr_long_df[lines_w_consistent_ids]
        # actually calculate the outliers
        expr_outlier_df = self.get_outliers(vcf_id_list)
        outs_per_id_file = re.sub('.txt', '_outliers_per_id_ALL',
                                  self.expr_outs_loc)
        plot_outs_per_id(expr_outlier_df, outs_per_id_file)
        outs_per_id_file = re.sub('.txt', '_outliers_per_id',
                                  self.expr_outs_loc)
        # determine which IDs have too many outliers (and remove these)
        if outlier_max:
            ids_to_keep = self.get_ids_w_low_out_ct(
                expr_outlier_df, outlier_max)
            lines_w_consistent_ids = self.expr_long_df.blinded_id.isin(
                ids_to_keep)
            if lines_w_consistent_ids.shape[0] == 0:
                raise RNASeqError("No IDs with less than {} outliers".format(
                    outlier_max))
            self.expr_long_df = self.expr_long_df[lines_w_consistent_ids]
            expr_outlier_df = self.get_outliers(ids_to_keep)
            plot_outs_per_id(expr_outlier_df, outs_per_id_file)
        # write `expr_outlier_df` to file
        print("Saving outlier status dataframe to", self.expr_outs_loc)
        expr_outlier_df.to_csv(self.expr_outs_loc, sep="\t", index=False)

    def get_outliers(self, ids_to_keep):
        """Calculate RNAseq outliers.

        Returns:
            `expr_outlier_df` (:obj:`DataFrame`): outliers per gene across
                all genes in long format

        Raises:
            :obj:`RNASeqError`: if there are no overlapping IDs between
                the RNAseq and WGS or if `distribution` is not valid

        Loads RNAseq in BED format (same format used for FastQTL),
            then identifies outliers.

        TODO:
            Clean/incorporate parallelization method and/or ask stackoverflow
                on options for vectorizing the `find_expr_outlier` function
            Confirm expr_cut_off is an acceptable value for normal/rank/custom
            Check that WGS IDs are consistent with outlier IDs

        """
        if self.distribution == "normal":
            self.identify_outliers_from_normal()
        elif self.distribution == "rank":
            self.identify_outliers_from_ranks(
                self.expr_long_df, self.least_extr_threshold)
        elif self.distribution == "custom":
            not_0_1 = ~self.expr_long_df.z_expr.isin([0, 1])
            if any(not_0_1):
                print(self.expr_long_df[not_0_1].head())
                raise RNASeqError("The values above were not 0 or 1")
            self.expr_long_df["expr_outlier"] = self.expr_long_df.z_expr == 1
            # set expr_outlier_neg as 0 for custom
            self.expr_long_df["expr_outlier_neg"] = 0
        else:
            raise RNASeqError("'{}' is not a valid outlier distribution".
                              format(self.distribution))
        if self.extrema:
            expr_outlier_df = self.find_most_extreme_expr_outlier()
        else:
            expr_outlier_df = self.expr_long_df
        expr_outlier_df.reset_index(inplace=True)
        return expr_outlier_df

    def identify_outliers_from_normal(self):
        """Identify outliers more extreme than a z-score threshold.

        TODO:
            All three lines raise a SettingWithCopyWarning when the column
            already exists in the dataframe. Unclear why

        """
        print("Calculating z-score outliers....")
        self.expr_long_df.loc[:, "z_abs"] = abs(self.expr_long_df.z_expr)
        self.expr_long_df.loc[:, "expr_outlier"] = (
            self.expr_long_df.z_abs > self.least_extr_threshold)
        self.expr_long_df.loc[:, "expr_outlier_neg"] = (
            (self.expr_long_df.z_expr < 0) &
            self.expr_long_df.expr_outlier)

    @staticmethod
    def identify_outliers_from_ranks(expr_long_df, least_extr_threshold):
        """Identify outliers based on those more extreme than percentile.

        Args
            `least_extr_threshold`: percentile cut-off for outliers

        """
        print("Calculating ranks...")
        expr_long_df = applyParallel(expr_long_df.groupby(
            'gene'), Outliers.calculate_ranks)
        print("Ranks calculated, identifying outliers")
        min_expr_cut_off = min(set(expr_long_df.expr_rank))
        if (least_extr_threshold <= min_expr_cut_off) or (
                least_extr_threshold >= 0.5):
            raise RNASeqError("The percentile cut-off specified ({}) is " +
                              "not between 0.5 and the minimum cut-off " +
                              "for this sample size, {}".format(
                                least_extr_threshold, min_expr_cut_off))
        hi_expr_cut_off = 1 - least_extr_threshold
        expr_long_df["expr_outlier_neg"] = (
            expr_long_df.expr_rank <= least_extr_threshold)
        expr_long_df["expr_outlier"] = (
            (expr_long_df.expr_rank >= hi_expr_cut_off) |
            expr_long_df.expr_outlier_neg)
        return expr_long_df

    @staticmethod
    def calculate_ranks(gene_group):
        """Calculate ranks for each gene.

        Args
            `gene_group`: expression for all IDs for a single gene

        Returns
            `gene_group`: with expr_rank which is the percentile

        """
        gene_group["expr_rank"] = gene_group["z_expr"].rank(method='average',
                                                            pct=True)
        return gene_group

    def test_normality(self):
        """Check if each gene has normal distribution.

        TODO:
            Options include QQ-plots, shapiro-wilk test and others.

        """
        return True

    def find_most_extreme_expr_outlier(self):
        """Loop over every gene in parallel, find the most extreme outlier.

        Returns:
            `expr_outlier_df` (:obj:`DataFrame`): outliers per gene across
                all genes in long format

        """
        print("Identifying most extreme outlier per gene...")
        expr_outlier_df = applyParallel(self.expr_long_df.groupby(
            'gene'), self.find_most_extreme_expr_outlier_per_gene)
        expr_outlier_df['expr_outlier_neg'] = (
            expr_outlier_df.expr_outlier_neg &
            expr_outlier_df.expr_outlier)
        return expr_outlier_df

    @staticmethod
    def find_most_extreme_expr_outlier_per_gene(gene_group):
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
