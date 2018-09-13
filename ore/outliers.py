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
from functools import partial

import pandas as pd
import matplotlib

from .genes import RNASeqError
from .plotting import plot_outs_per_id
from .utils import applyParallel

matplotlib.use('Agg')

import statsmodels.api as sm


"""statsmodels/compat/pandas.py:56:
FutureWarning: The pandas.core.datetools module is deprecated and will be
removed in a future version. Please use the pandas.tseries module instead.
"""


class Outliers(object):
    """Methods and attributes of outliers."""

    def __init__(self, pheno_loc, output_prefix, outlier_postfix,
                 extrema, distribution, threshold, cov, exclude_ids,
                 n_processes, logger):
        """Initialize outlier dataframe.

        Args:
            pheno_loc (:obj:`str`): gene expression (phenotype) location
            output_prefix (:obj:`str`): file prefix for outputs
            outlier_postfix (:obj:`str`): file ending for outlier files
            extrema (:obj:`boolean`): T/F for using most extreme outlier
            distribution (:obj:`str`): type of outlier distribution considered
            threshold (:obj:`list`): list of outlier cut-off thresholds
            n_processes (:obj:`int`): number of workers/cores to run at a time
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
        gene_expr_df = gene_expr_df.iloc[:, 3:]
        # logger.debug(gene_expr_df.head())
        self.n_processes = n_processes
        gene_expr_df.rename(columns={gene_expr_df.columns[0]: "gene"},
                            inplace=True)
        if exclude_ids:
            logger.debug("Expr DF size before excluding IDs:")
            logger.debug(gene_expr_df.shape)
            exclude_ids_in_df = [
                i for i in exclude_ids if i in gene_expr_df.columns]
            if exclude_ids_in_df:
                gene_expr_df.drop(exclude_ids_in_df, axis=1, inplace=True)
            logger.debug("Expr DF size AFTER excluding IDs:")
            logger.debug(gene_expr_df.shape)
        # if calculating covariates, re-normalize
        self.cov = cov
        """Re-calculating the z-score (not sure if appropriate)
        # if self.cov:
        gene_expr_df = self.recalculate_Zscore(gene_expr_df)
        # """
        # Convert gene expression data frame from wide to long:
        self.expr_long_df = pd.melt(
            gene_expr_df,
            id_vars='gene',  # gene_expr_df.columns.values[0],  # 'gene',
            value_vars=gene_expr_df.columns[1:].tolist(),
            var_name='blinded_id',
            value_name='z_expr')
        # logger.debug(self.expr_long_df.head())
        # logger.debug(self.expr_long_df.shape)
        # set the output file location
        self.expr_outs_loc = (output_prefix + "_outliers.txt")
        if outlier_postfix:
            self.expr_outs_loc = outlier_postfix
            # self.expr_outs_loc = (output_prefix + "_" + outlier_postfix)
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
        # logger.debug(self.least_extr_threshold)

    def prepare_outliers(self, outlier_max, vcf_id_list, logger):
        """Obtain gene expression outliers.

        Args:
            outlier_max (:obj:`int`): maximum number of outliers per ID
            vcf_id_list (:obj:`list`): list of Blinded IDs with WGS
            logger (:obj:`logging object`): Current logger

        Check if expression outlier file already exists. If it does not,
            use `get_outliers` to obtain outliers using only IDs
            with WGS. Once these outliers are obtained, plot a histogram
            of the number of outliers per ID. If there is a maximum
            number of outliers per ID specified, then remove IDs that
            cross this threshold and call outliers again. Save the
            outlier dataframe to `expr_outs_loc`

        """
        if os.path.exists(self.expr_outs_loc):
            print("Already made outlier file " + self.expr_outs_loc)
            return "Already made outlier file " + self.expr_outs_loc
        # only work with IDs in WGS that are also in the RNAseq
        lines_w_consistent_ids = self.expr_long_df.blinded_id.isin(vcf_id_list)
        if lines_w_consistent_ids.sum() == 0:
            raise RNASeqError("No overlapping IDs between RNAseq and VCF")
        self.expr_long_df = self.expr_long_df[lines_w_consistent_ids]
        if self.cov:
            print("before regression:")
            print(self.expr_long_df.head())
            self.expr_long_df = self.regress_out_covarates(self.cov)
            print("AFTER regression:")
            print(self.expr_long_df.head())
        # logger.debug(self.expr_long_df.head())
        # logger.debug(self.expr_long_df.shape)
        # actually calculate the outliers
        self.get_outliers(vcf_id_list)
        outs_per_id_file = re.sub('.txt', '_outliers_per_id_ALL',
                                  self.expr_outs_loc)
        plot_outs_per_id(self.expr_long_df, outs_per_id_file)
        outs_per_id_file = re.sub('.txt', '_outliers_per_id',
                                  self.expr_outs_loc)
        # determine which IDs have too many outliers (and remove these)
        if outlier_max:
            outs_per_id = self.expr_long_df[[
                'blinded_id', 'expr_outlier']].groupby('blinded_id').sum()
            while any(outs_per_id.expr_outlier >= outlier_max):
                ids_to_keep = self.get_ids_w_low_out_ct(
                    self.expr_long_df, outlier_max)
                lines_w_consistent_ids = self.expr_long_df.blinded_id.isin(
                    ids_to_keep)
                if lines_w_consistent_ids.shape[0] == 0:
                    raise RNASeqError("No IDs with <{} outliers".format(
                        outlier_max))
                self.expr_long_df = self.expr_long_df[lines_w_consistent_ids]
                self.get_outliers(ids_to_keep)
                plot_outs_per_id(self.expr_long_df, outs_per_id_file)
                outs_per_id = self.expr_long_df[[
                    'blinded_id', 'expr_outlier']].groupby('blinded_id').sum()
                # print(any(outs_per_id.expr_outlier >= outlier_max))
        self.remove_divergent_genes()
        # write `self.expr_long_df` to file
        print("Saving outlier status dataframe to", self.expr_outs_loc)
        self.expr_long_df.to_csv(self.expr_outs_loc, sep="\t", index=False)

    def get_outliers(self, ids_to_keep):
        """Calculate RNAseq outliers.

        Updates:
            `expr_long_df` (:obj:`DataFrame`): outliers per gene across
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
            self.identify_outliers_from_normal(ids_to_keep)
        elif self.distribution == "rank":
            self.expr_long_df = self.identify_outliers_from_ranks()
        elif self.distribution == "custom":
            not_0_1 = ~self.expr_long_df.z_expr.isin([0, 1])
            if any(not_0_1):
                print(self.expr_long_df[not_0_1].head())
                raise RNASeqError("The values above were not 0 or 1")
            self.expr_long_df["expr_outlier"] = self.expr_long_df.z_expr == 1
            # set expr_outlier_neg and expr_outlier_pos as 0 for custom
            self.expr_long_df["expr_outlier_neg"] = 0
            self.expr_long_df["expr_outlier_pos"] = 0
        else:
            raise RNASeqError("'{}' is not a valid outlier distribution".
                              format(self.distribution))
        if self.extrema:
            self.find_most_extreme_expr_outlier()

    def identify_outliers_from_normal(self, ids_to_keep):
        """Identify outliers more extreme than a z-score threshold.

        TODO:
            All three lines raise a SettingWithCopyWarning when the column
            already exists in the dataframe. Unclear why or if this is an issue

        """
        # print("(Re)calculating z-scores per gene...")
        print("Calculating z-score outliers....")
        self.expr_long_df = self.expr_long_df.assign(
            z_abs=abs(self.expr_long_df.z_expr))
        self.expr_long_df = self.expr_long_df.assign(
            expr_outlier=self.expr_long_df.z_abs > self.least_extr_threshold)
        self.expr_long_df = self.expr_long_df.assign(
            expr_outlier_neg=(self.expr_long_df.z_expr < 0) &
            self.expr_long_df.expr_outlier)
        self.expr_long_df = self.expr_long_df.assign(
            expr_outlier_pos=(self.expr_long_df.z_expr > 0) &
            self.expr_long_df.expr_outlier)
        # self.remove_divergent_genes(ids_to_keep)

    def remove_divergent_genes(self):
        """Remove genes where more than 5% of genes are outliers."""
        # self.expr_long_df.set_index(['gene', 'blinded_id'], inplace=True)
        # print(self.expr_long_df.index.get_level_values(
        #      'gene').unique())
        uniq_ids = self.expr_long_df.blinded_id.unique()
        print("Removing genes where more than 5% are outliers across " +
              str(len(uniq_ids)) + " samples.")
        if (self.distribution == "normal") and self.extrema:
            self.expr_long_df = self.expr_long_df.assign(
                expr_outlier_NOT_extrema=self.expr_long_df.z_abs >
                self.least_extr_threshold)
            outs_per_gene_ct = self.expr_long_df.groupby(
                'gene')['expr_outlier_NOT_extrema'].transform('sum')
            self.expr_long_df.drop(
                ['expr_outlier_NOT_extrema'], axis=1, inplace=True)
        elif self.distribution == "rank":
            # Temporary just to confirm using same genes across all comparisons
            return None
            self.expr_long_df = self.expr_long_df.assign(
                expr_outlier_NOT_rank=abs(self.expr_long_df.z_expr) > 2)
            outs_per_gene_ct = self.expr_long_df.groupby(
                'gene')['expr_outlier_NOT_rank'].transform('sum')
            self.expr_long_df.drop(
                ['expr_outlier_NOT_rank'], axis=1, inplace=True)
        else:
            outs_per_gene_ct = self.expr_long_df.groupby(
                'gene')['expr_outlier'].transform('sum')
        outs_per_gene_NOT_reasonable = (
            0.05*len(uniq_ids)) < outs_per_gene_ct
        # genes_to_rm = self.expr_long_df[
        #     outs_per_gene_NOT_reasonable].index.get_level_values(
        #     'gene').unique()
        genes_to_rm = self.expr_long_df[
            outs_per_gene_NOT_reasonable]['gene'].unique()
        print("More than 1/20 samples have outliers more more extreme " +
              "than Z={} for {} genes".format(
                  str(self.least_extr_threshold), str(len(genes_to_rm))))
        self.expr_long_df = self.expr_long_df[~outs_per_gene_NOT_reasonable]
        if self.expr_long_df.shape[0] == 0:
            raise RNASeqError("All genes have >1/20 samples as outliers")

    def identify_outliers_from_ranks(self):
        """Identify outliers based on those more extreme than percentile.

        Args
            `least_extr_threshold`: percentile cut-off for outliers

        """
        print("Calculating ranks...")
        expr_long_df = applyParallel(self.expr_long_df.groupby(
            'gene'), self.calculate_ranks,
            self.n_processes)
        print("Ranks calculated, identifying outliers")
        min_expr_cut_off = min(set(expr_long_df.expr_rank))
        if (self.least_extr_threshold <= min_expr_cut_off) or (
                self.least_extr_threshold >= 0.5):
            raise RNASeqError("The percentile cut-off specified ({}) is " +
                              "not between 0.5 and the minimum cut-off " +
                              "for this sample size, {}".format(
                                self.least_extr_threshold, min_expr_cut_off))
            # print(("The percentile cut-off specified ({}) is " +
            #        "not between 0.5 and the minimum cut-off " +
            #        "for this sample size, {}").format(
            #        self.least_extr_threshold, min_expr_cut_off))
            # self.least_extr_threshold = min_expr_cut_off
        hi_expr_cut_off = 1 - self.least_extr_threshold
        expr_long_df["expr_outlier_neg"] = (
            expr_long_df.expr_rank <= self.least_extr_threshold)
        expr_long_df["expr_outlier_pos"] = (
            expr_long_df.expr_rank >= hi_expr_cut_off)
        expr_long_df["expr_outlier"] = (
            expr_long_df.expr_outlier_pos |
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

    @staticmethod
    def recalculate_Zscore(expr_df):
        """Re-calculate the z-score for expression data."""
        expr_df.set_index(['gene'], inplace=True)
        expr_df_mean = expr_df.mean(axis=1)
        expr_df_std = expr_df.std(axis=1)
        expr_df = expr_df.sub(expr_df_mean, axis=0)
        expr_df = expr_df.div(expr_df_std, axis=0)
        expr_df.reset_index(inplace=True)
        return expr_df

    def regress_out_covarates(self, cov_loc):
        """Regress out covariates from the re-scaled expression matrix.

        Should we perform separate regressions for every gene or a single
            regression for all genes? Separate regressions because the
            parameters and variances will likely have wildly different
            estimates for different genes

        Source: https://stackoverflow.com/a/32102764

        """
        cov_df = pd.read_table(cov_loc, header=None)
        id_var = cov_df.iloc[0][0]
        cov_long = cov_df.set_index([0]).transpose().set_index(id_var)
        cov_long = cov_long.apply(pd.to_numeric, errors='ignore')
        calculate_residuals_per_gene_partial = partial(
            self.calculate_residuals_per_gene, cov_long=cov_long)
        expr_long_df_residuals = applyParallel(self.expr_long_df.groupby(
            'gene'), calculate_residuals_per_gene_partial,
            self.n_processes)
        return expr_long_df_residuals

    @staticmethod
    def calculate_residuals_per_gene(per_gene_df, cov_long):
        """Calculate the residuals per gene."""
        # make sure they have the same index before joining
        per_gene_df.set_index('blinded_id', inplace=True)
        current_gene = per_gene_df['gene'][0]
        del per_gene_df['gene']
        cov_long.index.name = per_gene_df.index.name
        per_gene_df.columns = ['gene_expr']
        # join covariates with expression
        per_gene_df = cov_long.join(per_gene_df, how='inner')
        # calculate residuals after regressing out covariates
        # sources: https://stackoverflow.com/a/32103366
        x_df = sm.add_constant(per_gene_df.iloc[:, :-1])
        model = sm.OLS(per_gene_df.gene_expr, x_df).fit()
        # return residuals added to the mean
        res_df = pd.DataFrame({'z_expr': model.resid + model.params.const})
        res_df.reset_index(inplace=True)
        res_df['gene'] = current_gene
        # return model.resid + model.params.const
        return res_df

    def test_normality(self):
        """Check if each gene has normal distribution.

        TODO:
            Options include QQ-plots, shapiro-wilk test and others.

        """
        return True

    def find_most_extreme_expr_outlier(self):
        """Loop over every gene in parallel, find the most extreme outlier.

        Updates attributes:
            `expr_long_df` (:obj:`DataFrame`): outliers per gene across
                all genes in long format

        """
        print("Identifying most extreme outlier per gene...")
        print(self.expr_long_df.head())
        print(self.expr_long_df.shape)
        self.expr_long_df = applyParallel(self.expr_long_df.groupby(
            'gene'), self.find_most_extreme_expr_outlier_per_gene,
            self.n_processes)
        self.expr_long_df['expr_outlier_neg'] = (
            self.expr_long_df.expr_outlier_neg &
            self.expr_long_df.expr_outlier)
        self.expr_long_df['expr_outlier_pos'] = (
            self.expr_long_df.expr_outlier_pos &
            self.expr_long_df.expr_outlier)

    @staticmethod
    def find_most_extreme_expr_outlier_per_gene(gene_group):
        """Label outliers in a group.

        Create a column for absolute value of expression z-score, determine
            which `blinded_id` has the maximum expression z-score to
            create `expr_outlier_status` column, then determine if any of these
            outliers have z-scores < 0 (i.e., are low or negative outliers).

        Args:
            `gene_group`: long-format expression dataframe for a gene

        """
        gene_group['expr_outlier'] = (
            (gene_group.z_abs == max(gene_group.z_abs)) &
            gene_group.expr_outlier)
        return gene_group

    @staticmethod
    def get_ids_w_low_out_ct(expr_outlier_df, outlier_max):
        """Identify and remove blinded_ids with a ton of outliers.

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
