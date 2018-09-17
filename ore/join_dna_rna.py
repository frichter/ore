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
import re

import pandas as pd

from .utils import filter_variant_class


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
                 annotations, contigs, logger):
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
        logger.info("Loading outliers...")
        self.load_outliers(expr_outs_loc)
        if annotations:
            anno_list = annotations
            anno_list = [re.sub("Conserved_TF_sites/", "Conserved_TF_", i)
                         for i in anno_list]
            anno_list = [re.sub("TfbsClustered_split/", "TfbsClust_", i)
                         for i in anno_list]
            anno_list = [re.sub("factorbookMotif/", "factorbookMotif_", i)
                         for i in anno_list]
            rep_w_blank = ".*/|.merged.sorted|.sorted|.bed$|.bed.gz$|.txt$"
            self.anno_list = [re.sub(rep_w_blank, "", i)
                              for i in anno_list]
        else:
            self.anno_list = None
        if os.path.exists(dna_rna_df_loc):
            logger.info("Already joined data, now loading from " +
                        dna_rna_df_loc)
            dtype_specs = {
                'dist_refgene': 'str', 'exon_func_refgene': 'str',
                'dist_ensgene': 'str', 'exon_func_ensgene': 'str'}
            self.df = pd.read_table(dna_rna_df_loc, dtype=dtype_specs)
            # self.df = pd.read_table(dna_rna_df_loc)
            self.df = filter_variant_class(self.df, variant_class, exon_class,
                                           refgene, ensgene)
            if variant_class:
                logger.info("Considering vars in these ENSEMBL categories " +
                            ", ".join(set(self.df.func_ensgene)))
            if exon_class:
                logger.info("Only vars in these ENSEMBL EXONIC categories " +
                            ", ".join(set(self.df.exon_func_ensgene)))
        else:
            logger.info("Loading variants...")
            self.load_vars(var_loc, contigs, variant_class, exon_class,
                           refgene, ensgene, max_tss_dist, logger)
            logger.info("joining outliers with variants...")
            # confirm there are overlapping IDs
            # logger.debug(self.var_df.head())
            # logger.debug(self.var_df.index[:10])
            # logger.debug(self.var_df.shape)
            # logger.debug(self.expr_outlier_df.head())
            # logger.debug(self.expr_outlier_df.index[:10])
            # logger.debug(self.expr_outlier_df.shape)
            dna_ids = self.var_df.index.levels[1]
            pheno_ids = self.expr_outlier_df.index.levels[1]
            overlapped_ids = dna_ids.isin(pheno_ids)
            if overlapped_ids.sum() == 0:
                raise JoinedDNAandRNAError("No overlapping IDs between" +
                                           "RNAseq and VCF")
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
        if self.anno_list:
            cols_to_keep.extend(self.anno_list)
        dtype_specs = {
            'dist_refgene': 'str', 'exon_func_refgene': 'str',
            'dist_ensgene': 'str', 'exon_func_ensgene': 'str'}
        for chrom in contigs:
            logger.info("Current chrom: " + chrom)  # "chr" +
            var_df_per_chrom = pd.read_table(
                var_loc % (chrom), dtype=dtype_specs)
            if var_df_per_chrom.empty:
                logger.info("Empty dataframe for chromosome " + chrom)
                raise JoinedDNAandRNAError("Empty dataframe for " + chrom)
            var_df_per_chrom.set_index(['gene', 'blinded_id'], inplace=True)
            var_df_per_chrom = var_df_per_chrom.loc[
                abs(var_df_per_chrom.tss_dist) <= max_tss_dist]
            var_df_per_chrom = filter_variant_class(
                var_df_per_chrom, variant_class, exon_class, refgene, ensgene)
            """# [18:118] [118:218] [218:-3]
            # last one is regions_enh_E013, total length is 371
            if len(cols_to_keep) == 14:
                cols_to_keep.extend(list(var_df_per_chrom)[18:118])
            logger.info(cols_to_keep)
            logger.info("Keeping {} columns".format(len(cols_to_keep)))
            # modification for summing accross annotations
            if 'any_gata4' in cols_to_keep:
                var_df_per_chrom = self.summarise_anno_cols(var_df_per_chrom)
            """
            var_df_per_chrom = var_df_per_chrom.reindex(columns=cols_to_keep)
            """Keep only lines where any of the annotations are 1.
            https://stackoverflow.com/a/34243246"""
            if self.anno_list:
                nzv = var_df_per_chrom[self.anno_list].any(axis=1)
                # print(nzv.sum())
                # print(var_df_per_chrom[nzv].shape)
                var_df_per_chrom = var_df_per_chrom[nzv]
            list_.append(var_df_per_chrom)
            # list_ is just pointing to DF so doesn't take any memory
            print(sys.getsizeof(var_df_per_chrom)/(1024**3), "Gb")
        logger.info("All contigs/chromosomes loaded")
        self.var_df = pd.concat(list_)
        logger.info(self.var_df.shape)
        if variant_class:
            logger.info("Considering variants in these ENSEMBL categories " +
                        ", ".join(set(self.var_df.func_ensgene)))
        if exon_class:
            logger.info("Only variants in these ENSEMBL EXONIC categories " +
                        ", ".join(set(self.var_df.exon_func_ensgene)))

    @staticmethod
    def summarise_anno_cols(df):
        """Keep a few prespecified summary columns for DeepHeart annos."""
        any_gata4 = [col for col in df.columns if 'ATA4' in col]  # ata4
        any_nkx25 = [col for col in df.columns if 'KX2' in col]  # kx2
        # any_tbx5 = [col for col in df.columns if 'TBX5' in col]  # bx5
        any_ep300 = [col for col in df.columns if 'EP300' in col]  # bx5
        any_polr2a = [col for col in df.columns if 'POLR2A' in col]  # bx5
        any_tbx = [col for col in df.columns if 'TBX' in col]  # bx5
        tf_search = re.compile('BRD4|EZH2|CTCF|SMARCA4')
        other_tf = [col for col in df.columns if tf_search.search(col)]  # bx5
        print(any_gata4)
        print(any_nkx25)
        # print(any_tbx5)
        print(any_ep300)
        print(any_polr2a)
        print(any_tbx)
        print(other_tf)
        all_tf = (any_gata4 + any_nkx25 + any_ep300 + any_polr2a +  # any_tbx5
                  any_tbx + other_tf)
        df['any_gata4'] = df[any_gata4].sum(axis=1) > 0
        df['any_nkx25'] = df[any_nkx25].sum(axis=1) > 0
        # df['any_tbx5'] = df[any_tbx5].sum(axis=1) > 0
        df['any_ep300'] = df[any_ep300].sum(axis=1) > 0
        df['any_polr2a'] = df[any_polr2a].sum(axis=1) > 0
        df['any_tbx'] = df[any_tbx].sum(axis=1) > 0
        df['other_tf'] = df[other_tf].sum(axis=1) > 0
        df['all_tf'] = df[all_tf].sum(axis=1) > 0
        # df['cvdc_enh_OR_prom'] = df[
        #     ['cvdc_enhancers_dickel',
        #      'cvdc_promoters.lineID']].sum(axis=1) > 0
        # explicitly convert to integers
        df.any_gata4 = df.any_gata4.astype(int)
        df.any_nkx25 = df.any_nkx25.astype(int)
        # df.any_tbx5 = df.any_tbx5.astype(int)
        df.any_ep300 = df.any_ep300.astype(int)
        df.any_polr2a = df.any_polr2a.astype(int)
        df.any_tbx = df.any_tbx.astype(int)
        df.other_tf = df.other_tf.astype(int)
        df.all_tf = df.all_tf.astype(int)
        # df.cvdc_enh_OR_prom = df.cvdc_enh_OR_prom.astype(int)
        return df

    def load_outliers(self, expr_outs_loc):
        """Load expression outlier dataframe.

        Attributes:
            expr_outlier_df (:obj:`DataFrame`): all gene-ID pairs with
                outliers labelled

        """
        self.expr_outlier_df = pd.read_table(expr_outs_loc)
        # self.expr_outlier_df = self.expr_outlier_df.iloc[:, 1:]
        self.expr_outlier_df.set_index(['gene', 'blinded_id'], inplace=True)

    def write_to_file(self, dna_rna_df_loc):
        """Write full joined DF to a file."""
        self.df.rename(columns={"var_id_freq": "intra_cohort_af",
                                "var_id_count": "intra_cohort_ac"},
                       inplace=True)
        cols_to_keep = ["blinded_id", "gene",
                        "gene_refgene", "func_refgene", "dist_refgene",
                        "exon_func_refgene", "gene_ensgene", "func_ensgene",
                        "dist_ensgene", "exon_func_ensgene",
                        "z_expr",
                        "expr_outlier", "expr_outlier_neg", "expr_outlier_pos",
                        "tss_dist", "var_id", "popmax_af", "intra_cohort_af",
                        "intra_cohort_ac", "VCF_af"]
        if self.anno_list:
            cols_to_keep.extend(self.anno_list)
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
