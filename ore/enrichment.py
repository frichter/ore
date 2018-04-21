#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Calculate enrichment.

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-02-01
:Copyright: 2018, Felix Richter
:License: CC BY-SA
"""

import itertools
import copy
from functools import partial
from multiprocessing import Pool, cpu_count

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact


class Enrich(object):
    """Methods and objects for enrichment.

    TODO:
        Separate the process of joining the final dataframe and
        calculating enrichment into different modules

    """

    def __init__(self, var_loc, expr_outs_loc, enrich_loc, rv_outlier_loc,
                 distribution, annovar_func, contigs):
        """Load and join variants and outliers.

        Args:
            `var_loc`: file location for the final set of variants that are
                being joined with expression (with string rep for chrom)
            `expr_outs_loc`: file location of the outliers
            `enrich_loc`: output file for enrichment calculations
            `rv_outlier_loc`: file location with the final set of
                outlier-variant pairs
            `distribution`: type of distribution being used for outliers
            `annovar_func`: annovar variant class to filter on (default None)
            `contigs`: chromosomes that are in the VCF

        Attributes:
            `joined_df`: variants and outliers in a single dataframe

        """
        self.var_loc = var_loc
        self.expr_outs_loc = expr_outs_loc
        self.enrich_loc = enrich_loc
        self.annovar_func = annovar_func
        if annovar_func:
            print("Keeping only {} variants".format(self.annovar_func))
        self.distribution = distribution
        self.rv_outlier_loc = rv_outlier_loc
        self.load_vars(contigs)
        self.load_outliers()
        print("joining outliers with variants...")
        self.joined_df = self.var_df.join(self.expr_outlier_df, how='inner')
        self.joined_df.reset_index(inplace=True)

    def load_vars(self, contigs):
        """Load and combine variant data.

        Attributes:
            `var_df`: all processed variants in a single dataframe

        """
        print("Loading variants...")
        self.var_df = pd.DataFrame()
        list_ = []
        cols_to_keep = ['popmax_af', 'var_id', 'tss_dist', 'annovar_func',
                        'var_id_freq']
        for chrom in contigs:
            # "wgs_pcgc_singletons_per_chrom/enh_var_hets_chr" + chrom + ".txt"
            print("chr" + chrom)
            var_df_per_chrom = pd.read_table(self.var_loc % ("chr" + chrom))
            # print("gene count:", len(var_df_per_chrom.gene.unique()))
            var_df_per_chrom.set_index(['gene', 'blinded_id'], inplace=True)
            # remove regions in repeats
            # var_df_per_chrom = var_df_per_chrom[var_df_per_chrom.rmsk == 0]
            if self.annovar_func:
                # vars_sub = var_df_per_chrom.annovar_func == self.annovar_func
                vars_sub = var_df_per_chrom.annovar_func.str.contains(
                    self.annovar_func)
                var_df_per_chrom = var_df_per_chrom[vars_sub]
            var_df_per_chrom = var_df_per_chrom[cols_to_keep]
            list_.append(var_df_per_chrom)
        self.var_df = pd.concat(list_)
        if self.annovar_func:
            print("Considering variants in the following refseq categories",
                  set(self.var_df.annovar_func))

    def load_outliers(self):
        """Load expression outlier dataframe.

        Attributes:
            `expr_outlier_df`: all gene-ID pairs with outliers labeled

        """
        print("Loading outliers...")
        self.expr_outlier_df = pd.read_table(self.expr_outs_loc)
        self.expr_outlier_df = self.expr_outlier_df.iloc[:, 1:]
        self.expr_outlier_df.set_index(['gene', 'blinded_id'], inplace=True)

    def loop_enrichment(self, n_processes, expr_cut_off_vec,
                        tss_cut_off_vec, af_cut_off_vec):
        """Loop over enrichment.

        Args:
            n_processes (:obj:`int`): number of processes to use

        """
        if isinstance(expr_cut_off_vec, float):
            expr_cut_off_vec = [expr_cut_off_vec]
        if isinstance(tss_cut_off_vec, float):
            tss_cut_off_vec = [tss_cut_off_vec]
        if isinstance(af_cut_off_vec, float):
            af_cut_off_vec = [af_cut_off_vec]
        cartesian_iter = itertools.product(expr_cut_off_vec,
                                           tss_cut_off_vec,
                                           af_cut_off_vec)
        # https://stackoverflow.com/questions/533905/get-the-cartesian-product-of-a-series-of-lists
        enrichment_per_tuple_partial = partial(
            self.enrichment_per_tuple)
        # run either multi-core or single core
        # if run_multi_core:
        print("Using {} cores, less than all {} cores".format(
              n_processes, cpu_count()))
        with Pool(n_processes) as p:
            out_line_list = p.map(enrichment_per_tuple_partial,
                                  cartesian_iter)
        # else:
        #     print("Running on a single process/core")
        #     out_line_list = []
        #     for cartesian_tuple in cartesian_iter:
        #         out_line_list.append(
        #             enrichment_per_tuple_partial(cartesian_tuple))
        self.write_enrichment_to_file(out_line_list)

    def enrichment_per_tuple(self, cut_off_tuple):
        """Calculate enrichment for each tuple.

        Deep copy the joined_df to avoid race conditions. Shallow copy is
            probably enough but deepcopy just to be safe

        """
        print("Calculating enrichment for", cut_off_tuple)
        enrich_df = copy.deepcopy(self.joined_df)
        max_intrapop_af = self.get_max_intra_pop_af(
            enrich_df, cut_off_tuple[2])
        enrich_df = self.identify_rows_to_keep(
            enrich_df, max_intrapop_af,
            self.distribution, cut_off_tuple)
        enrich_df = self.subset_deepcopy_df(enrich_df)
        # var_list = self.calculate_var_enrichment(enrich_df)
        gene_list = self.calculate_gene_enrichment(enrich_df)
        out_list = list(cut_off_tuple)
        # out_list.extend(var_list)
        out_list.extend(gene_list)
        print(out_list)
        return "\t".join([str(i) for i in out_list])

    def subset_deepcopy_df(self, enrich_df):
        """Subset a deep copy of the dataframe.

        Only keep genes with at least 1 rare variant.

        """
        enrich_df = enrich_df.loc[
            enrich_df.near_TSS & enrich_df.gene_has_out_w_vars]
        # confirm each gene has at least 1 rare variant
        # This line raises a SettingWithCopyWarning
        enrich_df.loc[:, 'gene_has_rare_vars'] = enrich_df.groupby(
            'gene')['rare_variant_status'].transform('sum') > 0
        enrich_df = enrich_df.loc[enrich_df.gene_has_rare_vars]
        return enrich_df

    @staticmethod
    def identify_rows_to_keep(joined_df, max_intrapop_af, distribution,
                              cut_off_tuple):
        """Reclassify rare variants and outliers and identify rows to keep.

        Keep only rows for genes with at least 1 outlier

        """
        expr_cut_off, tss_cut_off, af_cut_off = cut_off_tuple
        print("Parameters", expr_cut_off, tss_cut_off, af_cut_off)
        # classify as within x kb of TSS
        joined_df["near_TSS"] = abs(joined_df.tss_dist) < tss_cut_off
        # update outliers based on more extreme cut-offs
        if distribution == "normal":
            joined_df.expr_outlier = (
                joined_df.z_abs > expr_cut_off) & joined_df.expr_outlier
            joined_df.expr_outlier_neg = (joined_df.expr_outlier &
                                          joined_df.expr_outlier_neg)
        elif distribution == "rank":
            hi_expr_cut_off = 1 - expr_cut_off
            joined_df.loc[:, "expr_outlier_neg"] = (
                (joined_df.expr_rank <= expr_cut_off) &
                joined_df.expr_outlier_neg)
            joined_df.loc[:, "expr_outlier"] = (
                (joined_df.expr_rank >= hi_expr_cut_off) |
                joined_df.expr_outlier_neg) & joined_df.expr_outlier
        # else: raise error
        # only keep if there is any outlier in the gene
        joined_df.loc[:, 'gene_has_out_w_vars'] = joined_df.groupby(
            'gene')['expr_outlier'].transform('sum') > 0
        joined_df.loc[:, 'gene_has_NEG_out_w_vars'] = joined_df.groupby(
            'gene')['expr_outlier_neg'].transform('sum') > 0
        # classify as rare/common
        joined_df.loc[:, "rare_variant_status"] = (
            joined_df.popmax_af <= af_cut_off) & (
            joined_df.var_id_freq <= max_intrapop_af)
        # Only keep genes with at least 1 rare variant.
        # joined_df = joined_df[
        #     joined_df.near_TSS & joined_df.gene_has_out_w_vars]
        # # confirm each gene has at least 1 rare variant
        # joined_df['gene_has_rare_vars'] = joined_df.groupby(
        #     ['gene', 'tissue'])['rare_variant_status'].transform('sum') > 0
        # joined_df = joined_df[joined_df.gene_has_rare_vars]
        return joined_df

    @staticmethod
    def calculate_var_enrichment(enrich_df):
        """Calculate variant-centric enrichment."""
        out_tb = pd.crosstab(enrich_df.rare_variant_status,
                             enrich_df.expr_outlier)
        print(out_tb)
        neg_out_tb = pd.crosstab(enrich_df.rare_variant_status,
                                 enrich_df.expr_outlier_neg)
        (fet_or, fet_p), (fet_or_neg, fet_p_neg) = (
            fisher_exact(out_tb), fisher_exact(neg_out_tb))
        return [fet_or, fet_p, fet_or_neg, fet_p_neg]

    @staticmethod
    def calculate_gene_enrichment(enrich_df):
        """Calculate gene-centric enrichment.

        Counts each gene-ID pair once (i.e., removes extra variants per
            gene-ID pair with `drop_duplicates`)

        """
        enrich_df['gene_has_rare_var'] = enrich_df.groupby(
            ['gene', 'blinded_id'])['rare_variant_status'].transform('sum') > 0
        enrich_df = enrich_df.loc[:, [
            'gene', 'blinded_id', 'expr_outlier', 'expr_outlier_neg',
            'gene_has_NEG_out_w_vars', 'gene_has_rare_var']
            ].drop_duplicates(keep='first')
        out_tb = pd.crosstab(enrich_df.gene_has_rare_var,
                             enrich_df.expr_outlier)
        print(out_tb)
        out_list = Enrich.flatten_crosstab(out_tb)
        # out_list = out_tb.values.flatten().tolist()
        # while len(out_list) < 4:
        #     out_list.append("NA")
        # only keep subset of genes with low expression outliers
        enrich_df = enrich_df[enrich_df.gene_has_NEG_out_w_vars]
        neg_out_tb = pd.crosstab(enrich_df.gene_has_rare_var,
                                 enrich_df.expr_outlier_neg)
        print(neg_out_tb)
        neg_out_list = Enrich.flatten_crosstab(neg_out_tb)
        # neg_out_list = neg_out_tb.values.flatten().tolist()
        # while len(neg_out_list) < 4:
        #     neg_out_list.append("NA")
        out_list.extend(neg_out_list)
        try:
            fet_or, fet_p = fisher_exact(out_tb)
        except ValueError:
            fet_or, fet_p = ("NA", "NA")
        try:
            fet_or_neg, fet_p_neg = fisher_exact(neg_out_tb)
        except ValueError:
            fet_or_neg, fet_p_neg = ("NA", "NA")
        out_list.extend([fet_or, fet_p, fet_or_neg, fet_p_neg])
        return out_list

    def write_enrichment_to_file(self, out_line_list):
        """Write the enrichment results to a file."""
        with open(self.enrich_loc, 'w') as enrich_f:
            header_list = ["expr_cut_off", "tss_cut_off", "af_cut_off",
                           # "var_or", "var_p", "var_neg_or", "var_neg_p",
                           "not_rare_not_out", "not_rare_out",
                           "rare_not_out", "rare_out",
                           "not_rare_not_out_neg", "not_rare_out_neg",
                           "rare_not_out_neg", "rare_out_neg",
                           "gene_or", "gene_p", "gene_neg_or", "gene_neg_p"]
            enrich_f.write("\t".join(header_list) + "\n")
            for out_line in out_line_list:
                enrich_f.write(out_line + "\n")

    @staticmethod
    def flatten_crosstab(out_tb):
        """Flatten the crosstab output list."""
        out_list = out_tb.values.flatten().tolist()
        while len(out_list) < 4:
            if len(out_tb) > 0:
                if out_tb.index == np.array([True]):
                    out_list = [0] + out_list
                else:
                    out_list.append(0)
            else:
                out_list.append(0)
        return out_list

    @staticmethod
    def get_max_intra_pop_af(df_w_af, af_cut_off):
        """Obtain the maximum wihtin population AF to classify variants as rare.

        Args:
            `df_w_af` (:obj:`DataFrame`): a DF with a `var_id_freq` column
            `af_cut_off` (:obj:`float`): current allele frequency cut-off

        """
        if af_cut_off < min(df_w_af.var_id_freq):
            max_intra_pop_af = min(df_w_af.var_id_freq)
        else:
            max_intra_pop_af = af_cut_off
        return max_intra_pop_af

    def write_rvs_w_outs_to_file(self, out_cut_off, tss_cut_off, af_cut_off):
        """Write outliers with specified cut-offs to a file.

        Args:
            `out_cut_off`

        """
        cut_off_tuple = (out_cut_off, tss_cut_off, af_cut_off)
        enrich_df = copy.deepcopy(self.joined_df)
        max_intrapop_af = self.get_max_intra_pop_af(enrich_df, af_cut_off)
        print("Intra-population AF cut-off for rare:", max_intrapop_af)
        outlier_df = self.identify_rows_to_keep(
            enrich_df, max_intrapop_af,
            distribution=self.distribution, cut_off_tuple=cut_off_tuple)
        outlier_df = self.subset_deepcopy_df(outlier_df)
        # outlier_df.to_csv("test_all_joined.txt", index=False, sep="\t")
        # only keep outliers with rare variants
        outlier_df = outlier_df.loc[
            outlier_df.rare_variant_status & outlier_df.expr_outlier]
        cols_to_keep = ["blinded_id", "gene", "z_expr", "tss_dist",
                        "var_id", "popmax_af", "var_id_freq"]
        outlier_df = outlier_df[cols_to_keep]
        outlier_df.rename(columns={"var_id_freq": "intra_cohort_af"},
                          inplace=True)
        outlier_df.to_csv(self.rv_outlier_loc, index=False, sep="\t")

#
#
#
#
