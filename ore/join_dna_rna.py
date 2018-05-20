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


class JoinedDF(object):
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
            variant_class (:obj:`str`): annovar variant class to filter
                on (default None)
            exon_class (:obj:`str`): annovar EXON class to filter
                on (default None)
            contigs (:obj:`list`): chromosomes that are in the VCF
            logger (:obj:`logging object`): Current logger

        Attributes:
            `joined_df`: variants and outliers in a single dataframe

        """
        if os.path.exists(self.expr_outs_loc):
            logger.info("Already joined data")
            joined_df = pd.read_table(dna_rna_df_loc)
            return joined_df
        self.expr_outs_loc = expr_outs_loc
        logger.info("Loading variants...")
        self.load_vars(var_loc, contigs, variant_class, exon_class, refgene,
                       ensgene, logger)
        logger.info("Loading outliers...")
        self.load_outliers()
        logger.info("joining outliers with variants...")
        self.joined_df = self.var_df.join(self.expr_outlier_df, how='inner')
        logger.info(self.joined_df.head())
        logger.info(self.joined_df.shape)
        self.joined_df.reset_index(inplace=True)
        self.write_to_file(dna_rna_df_loc)
        return self.joined_df

    def load_vars(self, var_loc, contigs, variant_class, exon_class,
                  refgene, ensgene, logger):
        """Load and combine variant data.

        Attributes:
            `var_df`: all processed variants in a single dataframe

        """
        self.var_df = pd.DataFrame()
        list_ = []
        cols_to_keep = ['popmax_af', 'var_id', 'tss_dist',
                        'var_id_count', 'var_id_freq']
        dtype_specs = {
            'dist_refgene': 'str', 'exon_func_refgene': 'str',
            'dist_ensgene': 'str', 'exon_func_ensgene': 'str'}
        for chrom in contigs:
            logger.info("chr" + chrom)
            var_df_per_chrom = pd.read_table(
                var_loc % ("chr" + chrom), dtype=dtype_specs)
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

    def load_outliers(self):
        """Load expression outlier dataframe.

        Attributes:
            `expr_outlier_df`: all gene-ID pairs with outliers labeled

        """
        self.expr_outlier_df = pd.read_table(self.expr_outs_loc)
        self.expr_outlier_df = self.expr_outlier_df.iloc[:, 1:]
        self.expr_outlier_df.set_index(['gene', 'blinded_id'], inplace=True)

    def write_to_file(self, dna_rna_df_loc):
        """Write full joined DF to a file."""
        self.joined_df.rename(columns={"var_id_freq": "intra_cohort_af"},
                              inplace=True)
        self.joined_df.to_csv(dna_rna_df_loc, index=False, sep="\t")
