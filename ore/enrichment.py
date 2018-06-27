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
# from multiprocessing import Pool, cpu_count

from .enrich_utils import calculate_gene_enrichment
# from .enrich_utils import calculate_var_enrichment


class Enrich(object):
    """Methods and objects for enrichment.

    TODO:
        Separate the process of joining the final dataframe and
        calculating enrichment into different modules

    """

    def __init__(self, joined_df, enrich_loc, rv_outlier_loc, distribution):
        """Load and join variants and outliers.

        Args:
            joined_df (:obj:`DataFrame`): variants and outliers
            enrich_loc (:obj:`str`): output file for enrichment calculations
            rv_outlier_loc (:obj:`str`): file location with the final set of
                outlier-variant pairs
            distribution (:obj:`str`): type of distribution for outliers

        Attributes:
            `joined_df`: variants and outliers in a single dataframe

        """
        self.joined_df = joined_df
        self.enrich_loc = enrich_loc
        self.distribution = distribution
        self.rv_outlier_loc = rv_outlier_loc

    def loop_enrichment(self, n_processes, expr_cut_off_vec,
                        tss_cut_off_vec, af_cut_off_vec,
                        af_vcf, intracohort_rare_ac):
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
            self.enrichment_per_tuple,
            af_vcf=af_vcf, intracohort_rare_ac=intracohort_rare_ac)
        # print("Using {} cores, less than all {} cores".format(
        #       n_processes, cpu_count()))
        # with Pool(n_processes) as p:
        #     out_line_list = p.map(enrichment_per_tuple_partial,
        #                           cartesian_iter)
        out_line_list = []
        for cut_off_tuple in cartesian_iter:
            print(cut_off_tuple)
            out_list = enrichment_per_tuple_partial(cut_off_tuple)
            print(out_list)
            out_line_list.append(out_list)
        self.write_enrichment_to_file(out_line_list)

    def enrichment_per_tuple(self, cut_off_tuple, af_vcf=False,
                             intracohort_rare_ac=None):
        """Calculate enrichment for each tuple.

        Deep copy the joined_df to avoid race conditions. Shallow copy is
            probably enough but deepcopy just to be safe

        """
        print("Calculating enrichment for", cut_off_tuple)
        enrich_df = copy.deepcopy(self.joined_df)
        # replace af_cut_off with intra-cohort minimum if former is
        # smaller than latter
        print(cut_off_tuple)
        max_intrapop_af = self.get_max_intra_pop_af(
            enrich_df, cut_off_tuple[2])
        max_vcf_af = self.get_max_vcf_af(enrich_df, cut_off_tuple[2])
        enrich_df = self.identify_rows_to_keep(
            enrich_df, max_intrapop_af, max_vcf_af,
            af_vcf, intracohort_rare_ac,
            self.distribution, cut_off_tuple)
        # var_list = calculate_var_enrichment(enrich_df)
        gene_list = calculate_gene_enrichment(enrich_df)
        out_list = list(cut_off_tuple)
        # out_list.extend(var_list)
        out_list.extend(gene_list)
        return "\t".join([str(i) for i in out_list])

    @staticmethod
    def identify_rows_to_keep(joined_df, max_intrapop_af, max_vcf_af,
                              af_vcf, intracohort_rare_ac,
                              distribution, cut_off_tuple):
        """Reclassify rare variants and outliers and identify rows to keep.

        Keep only rows for genes with at least 1 outlier

        """
        expr_cut_off, tss_cut_off, af_cut_off = cut_off_tuple
        print("Parameters", expr_cut_off, tss_cut_off, af_cut_off)
        # classify as within x kb of TSS
        joined_df["near_TSS"] = abs(joined_df.tss_dist) <= tss_cut_off
        joined_df = joined_df.loc[joined_df.near_TSS]
        # print("filtered by TSS:", joined_df.shape)
        # update outliers based on more extreme cut-offs
        if distribution == "normal":
            joined_df.expr_outlier = (
                joined_df.z_abs >= expr_cut_off) & joined_df.expr_outlier
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
        # expression outliers should be preset as 0 and 1 for custom
        # i.e., there should be no parameter iteration in this dimension
        # else: raise error
        # only keep if there is any outlier in the gene
        joined_df.loc[:, 'gene_has_out_w_vars'] = joined_df.groupby(
            'gene')['expr_outlier'].transform('sum') > 0
        joined_df.loc[:, 'gene_has_NEG_out_w_vars'] = joined_df.groupby(
            'gene')['expr_outlier_neg'].transform('sum') > 0
        joined_df = joined_df.loc[joined_df.gene_has_out_w_vars]
        # classify as rare/common
        rare_variant_status = joined_df.popmax_af <= af_cut_off
        if af_vcf:
            rare_variant_status = rare_variant_status & (
                joined_df.VCF_af <= max_vcf_af)
        if intracohort_rare_ac:
            rare_variant_status = rare_variant_status & (
                joined_df.intra_cohort_ac <= intracohort_rare_ac)
        else:
            rare_variant_status = rare_variant_status & (
                joined_df.intra_cohort_af <= max_intrapop_af)
        joined_df.loc[:, "rare_variant_status"] = rare_variant_status
        return joined_df

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
                if not out_line.startswith("NA_line"):
                    enrich_f.write(out_line + "\n")
                else:
                    print("Skipping", out_line)

    @staticmethod
    def get_max_intra_pop_af(df_w_af, af_cut_off):
        """Obtain the maximum wihtin population AF to classify variants as rare.

        Args:
            `df_w_af` (:obj:`DataFrame`): a DF with a `intra_cohort_af` column
            `af_cut_off` (:obj:`float`): current allele frequency cut-off

        """
        if af_cut_off < min(df_w_af.intra_cohort_af):
            max_intra_pop_af = min(df_w_af.intra_cohort_af)
        else:
            max_intra_pop_af = af_cut_off
        return max_intra_pop_af

    @staticmethod
    def get_max_vcf_af(df_w_vcf_af, af_cut_off):
        """Obtain the maximum VCF AF to classify variants as rare.

        Args:
            `df_w_vcf_af` (:obj:`DataFrame`): a DF with a `VCF_af` column
            `af_cut_off` (:obj:`float`): current allele frequency cut-off

        """
        if af_cut_off < min(df_w_vcf_af.VCF_af):
            max_vcf_af = min(df_w_vcf_af.VCF_af)
        else:
            max_vcf_af = af_cut_off
        return max_vcf_af

    def write_rvs_w_outs_to_file(self, out_cut_off, tss_cut_off, af_cut_off,
                                 af_vcf, intracohort_rare_ac):
        """Write outliers with specified cut-offs to a file.

        Args:
            `out_cut_off`

        """
        cut_off_tuple = (out_cut_off, tss_cut_off, af_cut_off)
        enrich_df = copy.deepcopy(self.joined_df)
        max_intrapop_af = self.get_max_intra_pop_af(enrich_df, af_cut_off)
        max_vcf_af = self.get_max_vcf_af(enrich_df, af_cut_off)
        outlier_df = self.identify_rows_to_keep(
            enrich_df, max_intrapop_af, max_vcf_af,
            af_vcf, intracohort_rare_ac,
            distribution=self.distribution, cut_off_tuple=cut_off_tuple)
        # only keep outliers with rare variants
        outlier_df = outlier_df.loc[
            outlier_df.rare_variant_status & outlier_df.expr_outlier]
        cols_to_keep = ["blinded_id", "gene", "z_expr", "tss_dist",
                        "var_id", "popmax_af", "intra_cohort_af",
                        "intra_cohort_ac", "VCF_af"]
        outlier_df = outlier_df[cols_to_keep]
        outlier_df.to_csv(self.rv_outlier_loc, index=False, sep="\t",
                          float_format='%g')

#
#
#
#
