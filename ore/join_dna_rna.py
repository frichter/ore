#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Join outlier (RNA) and variant (DNA) info.

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-05-20
:Copyright: 2018, Felix Richter
:License: CC BY-SA
"""


import os

import pandas as pd


class JoinedVarExpr(object):
    """Methods for creating the final DF."""

    def __init__(self, var_loc, expr_outs_loc, dna_rna_df_loc,
                 variant_class, exon_class, refgene, ensgene,
                 contigs, logger):
        """Load and join variants and outliers.

        Args:
            var_loc (:obj:`str`): file location for the final set of variants
                that are being joined with expression (contains a
                string rep for chrom)
            expr_outs_loc (:obj:`str`): file location of the outliers
            dna_rna_df_loc (:obj:`str`): fully joined dataframe output file
            variant_class (:obj:`str`): annovar variant class to filter
                on (default None)
            exon_class (:obj:`str`): annovar EXON class to filter
                on (default None)
            contigs (:obj:`list`): chromosomes that are in the VCF
            logger (:obj:`logging object`): Current logger

        Attributes:
            df (:obj:`DataFrame`): variants and outliers in a single dataframe

        """
        if os.path.exists(dna_rna_df_loc):
            logger.info("Already joined data")
            self.df = pd.read_table(dna_rna_df_loc)
        else:
            logger.info("Loading variants...")
            self.load_vars(var_loc, contigs, variant_class, exon_class,
                           refgene, ensgene, logger)
            logger.info("Loading outliers...")
            self.load_outliers(expr_outs_loc)
            logger.info("joining outliers with variants...")
            logger.debug(self.var_df.head())
            logger.debug(self.var_df.index[:10])
            logger.debug(self.var_df.shape)
            logger.debug(self.expr_outlier_df.head())
            logger.debug(self.expr_outlier_df.index[:10])
            logger.debug(self.expr_outlier_df.shape)
            self.df = self.var_df.join(self.expr_outlier_df, how='inner')
            self.df.reset_index(inplace=True)
            self.write_to_file(dna_rna_df_loc)

    def load_vars(self, var_loc, contigs, variant_class, exon_class,
                  refgene, ensgene, logger):
        """Load and combine variant data.

        Attributes:
            var_df (:obj:`DataFrame`): all variants in a single dataframe

        """
        self.var_df = pd.DataFrame()
        list_ = []
        cols_to_keep = ['popmax_af', 'var_id', 'tss_dist', 'VCF_af',
                        'var_id_count', 'var_id_freq']
        dtype_specs = {
            'dist_refgene': 'str', 'exon_func_refgene': 'str',
            'dist_ensgene': 'str', 'exon_func_ensgene': 'str'}
        for chrom in contigs:
            logger.info("Current chrom: " + chrom)  # "chr" +
            var_df_per_chrom = pd.read_table(
                var_loc % (chrom), dtype=dtype_specs)
            var_df_per_chrom.set_index(['gene', 'blinded_id'], inplace=True)
            if variant_class:
                var_df_per_chrom = self.filter_refgene_ensgene(
                    var_df_per_chrom, variant_class, refgene, ensgene)
                cols_to_keep.extend(['func_refgene', 'func_ensgene'])
            if exon_class:
                var_df_per_chrom = self.filter_refgene_ensgene_exon(
                    var_df_per_chrom, exon_class, refgene, ensgene)
                cols_to_keep.extend(['exon_func_refgene', 'exon_func_ensgene'])
            var_df_per_chrom = var_df_per_chrom[cols_to_keep]
            list_.append(var_df_per_chrom)
        logger.info("All contigs/chromosomes loaded")
        self.var_df = pd.concat(list_)
        logger.info(self.var_df.shape)
        if variant_class:
            logger.info("Considering variants in the following categories",
                        set(self.var_df.func_refgene))
        if exon_class:
            logger.info("Only variants in the following EXONIC categories",
                        set(self.var_df.exon_func_refgene))

    @staticmethod
    def filter_refgene_ensgene(var_df_per_chrom, variant_class,
                               refgene, ensgene):
        """Filter for a refgene function, ensembl function or both."""
        if refgene:
            vars_refgene = var_df_per_chrom.func_refgene.str.contains(
                variant_class)
            var_df_per_chrom = var_df_per_chrom[vars_refgene]
        if ensgene:
            vars_ensgene = var_df_per_chrom.func_ensgene.str.contains(
                variant_class)
            var_df_per_chrom = var_df_per_chrom[vars_ensgene]
        return var_df_per_chrom

    @staticmethod
    def filter_refgene_ensgene_exon(var_df_per_chrom, exon_class,
                                    refgene, ensgene):
        """Filter for a refgene function, ensembl function or both."""
        if refgene:
            vars_refgene = var_df_per_chrom.exon_func_refgene.str.contains(
                exon_class)
            var_df_per_chrom = var_df_per_chrom[vars_refgene]
        if ensgene:
            vars_ensgene = var_df_per_chrom.exon_func_ensgene.str.contains(
                exon_class)
            var_df_per_chrom = var_df_per_chrom[vars_ensgene]
        return var_df_per_chrom

    def load_outliers(self, expr_outs_loc):
        """Load expression outlier dataframe.

        Attributes:
            expr_outlier_df (:obj:`DataFrame`): all gene-ID pairs with
                outliers labelled

        """
        self.expr_outlier_df = pd.read_table(expr_outs_loc)
        self.expr_outlier_df = self.expr_outlier_df.iloc[:, 1:]
        self.expr_outlier_df.set_index(['gene', 'blinded_id'], inplace=True)

    def write_to_file(self, dna_rna_df_loc):
        """Write full joined DF to a file."""
        self.df.rename(columns={"var_id_freq": "intra_cohort_af",
                                "var_id_count": "intra_cohort_ac"},
                       inplace=True)
        print(self.df.head())
        print(self.df.shape)
        cols_to_keep = ["blinded_id", "gene", "z_expr", "z_abs",
                        "expr_outlier", "expr_outlier_neg", "tss_dist",
                        "var_id", "popmax_af", "intra_cohort_af",
                        "intra_cohort_ac", "VCF_af"]
        out_df = self.df[cols_to_keep]
        out_df.to_csv(dna_rna_df_loc, index=False, sep="\t",
                      float_format='%g')

#
#
