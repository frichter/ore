#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Join outlier (RNA) and variant (DNA) info.

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-05-20
:Copyright: 2018, Felix Richter
:License: CC BY-SA
"""


import os
import sys

import pandas as pd


class Error(Exception):
    """Base class for exceptions in this module."""

    pass


class JoinedDNAandRNAError(Error):
    """Exception raised for errors in this module.

    Attributes:
        message -- explanation of the error

    """

    def __init__(self, message):
        """Assign error explanation to object."""
        self.message = message


class JoinedVarExpr(object):
    """Methods for creating the final DF."""

    def __init__(self, var_loc, expr_outs_loc, dna_rna_df_loc,
                 variant_class, exon_class, refgene, ensgene, max_tss_dist,
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
            max_tss_dist (:obj:`int`): maximum distance from TSS
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
                           refgene, ensgene, max_tss_dist, logger)
            logger.info("Loading outliers...")
            self.load_outliers(expr_outs_loc)
            logger.info("joining outliers with variants...")
            # logger.debug(self.var_df.head())
            # logger.debug(self.var_df.index[:10])
            # logger.debug(self.var_df.shape)
            # logger.debug(self.expr_outlier_df.head())
            # logger.debug(self.expr_outlier_df.index[:10])
            # logger.debug(self.expr_outlier_df.shape)
            self.df = self.var_df.join(self.expr_outlier_df, how='inner')
            self.df.reset_index(inplace=True)
            self.write_to_file(dna_rna_df_loc)

    def load_vars(self, var_loc, contigs, variant_class, exon_class,
                  refgene, ensgene, max_tss_dist, logger):
        """Load and combine variant data.

        Attributes:
            var_df (:obj:`DataFrame`): all variants in a single dataframe

        """
        self.var_df = pd.DataFrame()
        list_ = []
        # cols_to_keep does not include 'gene' and 'blinded_id'
        # because these are set as the indices
        cols_to_keep = ['var_id', 'tss_dist',
                        'gene_refgene', 'func_refgene', 'dist_refgene',
                        'exon_func_refgene',
                        'gene_ensgene', 'func_ensgene', 'dist_ensgene',
                        'exon_func_ensgene',
                        'popmax_af', 'VCF_af', 'var_id_count', 'var_id_freq']
        dtype_specs = {
            'dist_refgene': 'str', 'exon_func_refgene': 'str',
            'dist_ensgene': 'str', 'exon_func_ensgene': 'str'}
        for chrom in contigs:
            logger.info("Current chrom: " + chrom)  # "chr" +
            var_df_per_chrom = pd.read_table(
                var_loc % (chrom), dtype=dtype_specs)
            var_df_per_chrom.set_index(['gene', 'blinded_id'], inplace=True)
            var_df_per_chrom = var_df_per_chrom.loc[
                abs(var_df_per_chrom.tss_dist) <= max_tss_dist]
            if variant_class:
                var_df_per_chrom = self.filter_refgene_ensgene(
                    var_df_per_chrom, variant_class, refgene, ensgene)
                # cols_to_keep.extend(['func_refgene', 'func_ensgene'])
                if variant_class.startswith("exonic") and exon_class:
                    var_df_per_chrom = self.filter_refgene_ensgene_exon(
                        var_df_per_chrom, exon_class, refgene, ensgene)
                    # cols_to_keep.extend(['exon_func_refgene',
                    #                      'exon_func_ensgene'])
            var_df_per_chrom = var_df_per_chrom[cols_to_keep]
            list_.append(var_df_per_chrom)
            print(sys.getsizeof(var_df_per_chrom)/(1024**3), "Gb")
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
                variant_class, regex=True)
            var_df_per_chrom = var_df_per_chrom[vars_refgene]
        if ensgene:
            vars_ensgene = var_df_per_chrom.func_ensgene.str.contains(
                variant_class, regex=True)
            var_df_per_chrom = var_df_per_chrom[vars_ensgene]
        return var_df_per_chrom

    @staticmethod
    def filter_refgene_ensgene_exon(var_df_per_chrom, exon_class,
                                    refgene, ensgene):
        """Filter for a refgene function, ensembl function or both.

        Args:
            var_df_per_chrom (:obj:`DataFrame`): all variants in a chromosome
            variant_class (:obj:`str`): annovar variant class to filter
                on (default None)
            exon_class (:obj:`str`): annovar EXON class to filter
                on (default None)
            refgene (:obj:`boolean`): if used RefSeq to define variant classes
            ensgene (:obj:`boolean`): using ENSEMBL to define variant classes

        Returns:
            var_df_per_chrom (:obj:`DataFrame`): only variants in the
                desired `exon_class`

        Description:
            First prepends a ^ so that only the highest impact `exon_class`
            is considered as the de-facto class for filtering.

        """
        exon_class = "^" + exon_class
        print(exon_class)
        if refgene:
            vars_refgene = var_df_per_chrom.exon_func_refgene.str.contains(
                exon_class, regex=True)
            print(var_df_per_chrom.head)
            print(var_df_per_chrom.shape)
            var_df_per_chrom = var_df_per_chrom[vars_refgene]
            print(var_df_per_chrom.head)
            print(var_df_per_chrom.shape)
        if ensgene:
            vars_ensgene = var_df_per_chrom.exon_func_ensgene.str.contains(
                exon_class, regex=True)
            var_df_per_chrom = var_df_per_chrom[vars_ensgene]
            print(var_df_per_chrom.head)
            print(var_df_per_chrom.shape)
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
        cols_to_keep = ["blinded_id", "gene", "z_expr",
                        "expr_outlier", "expr_outlier_neg", "tss_dist",
                        "var_id", "popmax_af", "intra_cohort_af",
                        "intra_cohort_ac", "VCF_af"]
        missing_cols = [i for i in cols_to_keep if i not in self.df.columns]
        if any(missing_cols):
            raise JoinedDNAandRNAError("The following essential columns are " +
                                       "missing: " + ", ".join(missing_cols))
        # keeping all columns for now
        # if "z_abs" in self.df.columns:
        #     cols_to_keep.insert(3, "z_abs")
        out_df = self.df  # [cols_to_keep]
        out_df.to_csv(dna_rna_df_loc, index=False, sep="\t",
                      float_format='%g')

#
#
