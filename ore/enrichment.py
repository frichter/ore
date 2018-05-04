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
                 distribution, variant_class, refgene, ensgene, contigs):
        """Load and join variants and outliers.

        Args:
            `var_loc`: file location for the final set of variants that are
                being joined with expression (with string rep for chrom)
            `expr_outs_loc`: file location of the outliers
            `enrich_loc`: output file for enrichment calculations
            `rv_outlier_loc`: file location with the final set of
                outlier-variant pairs
            `distribution`: type of distribution being used for outliers
            `variant_class`: annovar variant class to filter on (default None)
            `contigs`: chromosomes that are in the VCF

        Attributes:
            `joined_df`: variants and outliers in a single dataframe

        """
        self.var_loc = var_loc
        self.expr_outs_loc = expr_outs_loc
        self.enrich_loc = enrich_loc
        self.distribution = distribution
        self.rv_outlier_loc = rv_outlier_loc
        self.load_vars(contigs, variant_class, refgene, ensgene)
        self.load_outliers()
        print("joining outliers with variants...")
        self.joined_df = self.var_df.join(self.expr_outlier_df, how='inner')
        self.joined_df.reset_index(inplace=True)

    def load_vars(self, contigs, variant_class, refgene, ensgene):
        """Load and combine variant data.

        Attributes:
            `var_df`: all processed variants in a single dataframe

        """
        print("Loading variants...")
        self.var_df = pd.DataFrame()
        list_ = []
        cols_to_keep = ['popmax_af', 'var_id', 'tss_dist', 'func_refgene',
                        'exon_func_refgene', 'func_ensgene',
                        'exon_func_ensgene', 'var_id_count', 'var_id_freq']
        cols_to_keep.extend(['nkx2.5.mm9.hg19', 'regions_enh_E013'])
        cols_to_keep.extend(
            ['any_gata4', 'any_nkx25', 'any_tbx5', 'all_tf',
             'cvdc_enh_OR_prom'])
        dtype_specs = {
            'dist_refgene': 'str', 'exon_func_refgene': 'str',
            'dist_ensgene': 'str', 'exon_func_ensgene': 'str'}
        for chrom in contigs:  # ["21", "22"]:  #
            # "wgs_pcgc_singletons_per_chrom/enh_var_hets_chr" + chrom + ".txt"
            print("chr" + chrom)
            var_df_per_chrom = pd.read_table(
                self.var_loc % ("chr" + chrom), dtype=dtype_specs)
            var_df_per_chrom.set_index(['gene', 'blinded_id'], inplace=True)
            # remove regions in repeats
            if variant_class:
                var_df_per_chrom = self.filter_refgene_ensgene(
                    var_df_per_chrom, variant_class, refgene, ensgene)
            # [18:118] [118:218] [218:-3]
            # last one is regions_enh_E013, total length is 371
            # if len(cols_to_keep) == 9:
            #     cols_to_keep.extend(list(var_df_per_chrom)[18+325:-3])
            # modification for summing accross annotations
            var_df_per_chrom = self.summarise_anno_cols(var_df_per_chrom)
            var_df_per_chrom = var_df_per_chrom[cols_to_keep]
            print(var_df_per_chrom.head(2))
            print(var_df_per_chrom.all_tf.sum())
            print(var_df_per_chrom.any_nkx25.sum())
            list_.append(var_df_per_chrom)
        print("All contigs/chromosomes loaded")
        self.var_df = pd.concat(list_)
        print(self.var_df.shape)
        if variant_class:
            print("Considering variants in the following categories",
                  set(self.var_df.func_refgene))
        # if exon_class:
        #     print("Considering variants in the following Exonic categories",
        #           set(self.var_df.exon_func_refgene))

    @staticmethod
    def summarise_anno_cols(df):
        """Keep a few prespecified summary columns."""
        any_gata4 = [col for col in df.columns if 'ata4' in col]
        any_nkx25 = [col for col in df.columns if 'kx2' in col]
        any_tbx5 = [col for col in df.columns if 'bx5' in col]
        all_tf = any_gata4 + any_nkx25 + any_tbx5
        print(all_tf)
        df['any_gata4'] = df[any_gata4].sum(axis=1) > 0
        df['any_nkx25'] = df[any_nkx25].sum(axis=1) > 0
        df['any_tbx5'] = df[any_tbx5].sum(axis=1) > 0
        df['all_tf'] = df[all_tf].sum(axis=1) > 0
        df['cvdc_enh_OR_prom'] = df[
            ['cvdc_enhancers_dickel', 'cvdc_promoters.lineID']].sum(axis=1) > 0
        df.any_gata4 = df.any_gata4.astype(int)
        df.any_nkx25 = df.any_nkx25.astype(int)
        df.any_tbx5 = df.any_tbx5.astype(int)
        df.all_tf = df.all_tf.astype(int)
        df.cvdc_enh_OR_prom = df.cvdc_enh_OR_prom.astype(int)
        return df

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
        print("Anno column final index:", self.joined_df.shape[1] - 6)
        anno_list = [i for i in range(11, self.joined_df.shape[1] - 6)]
        if isinstance(expr_cut_off_vec, float):
            expr_cut_off_vec = [expr_cut_off_vec]
        if isinstance(tss_cut_off_vec, float):
            tss_cut_off_vec = [tss_cut_off_vec]
        if isinstance(af_cut_off_vec, float):
            af_cut_off_vec = [af_cut_off_vec]
        cartesian_iter = itertools.product(expr_cut_off_vec,
                                           tss_cut_off_vec,
                                           af_cut_off_vec,
                                           anno_list)
        # https://stackoverflow.com/questions/533905/get-the-cartesian-product-of-a-series-of-lists
        enrichment_per_tuple_partial = partial(
            self.enrichment_per_tuple)
        # print("Using {} cores, less than all {} cores".format(
        #       n_processes, cpu_count()))
        # with Pool(n_processes) as p:
        #     out_line_list = p.map(enrichment_per_tuple_partial,
        #                           cartesian_iter)
        out_line_list = []
        for cut_off_tuple in cartesian_iter:
            out_list = enrichment_per_tuple_partial(cut_off_tuple)
            print(out_list)
            out_line_list.append(out_list)
        self.write_enrichment_to_file(out_line_list)

    def enrichment_per_tuple(self, cut_off_tuple):
        """Calculate enrichment for each tuple.

        Deep copy the joined_df to avoid race conditions. Shallow copy is
            probably enough but deepcopy just to be safe

        """
        print("Calculating enrichment for", cut_off_tuple)
        enrich_df = copy.deepcopy(self.joined_df)
        current_anno = list(enrich_df)[cut_off_tuple[3]]
        print("current column:", current_anno)
        in_anno = enrich_df.loc[:, current_anno] == 1
        # keep only a specific annotation
        enrich_df = enrich_df.loc[in_anno]
        # remove annotation column index number from tuple
        cut_off_tuple = tuple(list(cut_off_tuple)[:-1])
        if enrich_df.shape[0] == 0:
            return "NA_line: no overlaps with " + current_anno
        # replace af_cut_off with intra-cohort minimum if former is
        # smaller than latter
        max_intrapop_af = self.get_max_intra_pop_af(
            enrich_df, cut_off_tuple[2])
        enrich_df = self.identify_rows_to_keep(
            enrich_df, max_intrapop_af,
            self.distribution, cut_off_tuple)
        # var_list = self.calculate_var_enrichment(enrich_df)
        gene_list = self.calculate_gene_enrichment(enrich_df)
        out_list = list(cut_off_tuple)
        # out_list.extend(var_list)
        out_list.extend(gene_list)
        return "\t".join([str(i) for i in out_list]) + "\t" + current_anno

    @staticmethod
    def identify_rows_to_keep(joined_df, max_intrapop_af, distribution,
                              cut_off_tuple):
        """Reclassify rare variants and outliers and identify rows to keep.

        Keep only rows for genes with at least 1 outlier

        """
        expr_cut_off, tss_cut_off, af_cut_off = cut_off_tuple
        print("Parameters", expr_cut_off, tss_cut_off, af_cut_off)
        # classify as within x kb of TSS
        joined_df["near_TSS"] = abs(joined_df.tss_dist) <= tss_cut_off
        joined_df = joined_df.loc[joined_df.near_TSS]
        print("filtered by TSS:", joined_df.shape)
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
        # else: raise error
        # only keep if there is any outlier in the gene
        joined_df.loc[:, 'gene_has_out_w_vars'] = joined_df.groupby(
            'gene')['expr_outlier'].transform('sum') > 0
        joined_df.loc[:, 'gene_has_NEG_out_w_vars'] = joined_df.groupby(
            'gene')['expr_outlier_neg'].transform('sum') > 0
        joined_df = joined_df.loc[joined_df.gene_has_out_w_vars]
        print("filtered by if gene has outlier:", joined_df.shape)
        # classify as rare/common
        joined_df.loc[:, "rare_variant_status"] = (
            joined_df.popmax_af <= af_cut_off) & (
            # joined_df.var_id_freq <= max_intrapop_af)
            joined_df.var_id_count <= 5)
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
                           "gene_or", "gene_p", "gene_neg_or", "gene_neg_p",
                           "annotaion"]
            enrich_f.write("\t".join(header_list) + "\n")
            for out_line in out_line_list:
                if not out_line.startswith("NA_line"):
                    enrich_f.write(out_line + "\n")
                else:
                    print("Skipping", out_line)

    @staticmethod
    def flatten_crosstab(out_tb):
        """Flatten the crosstab output list."""
        out_list = out_tb.values.flatten().tolist()
        while len(out_list) < 4:
            # if there's only 1 category in gene_has_rare_var...
            if len(out_tb) == 1:
                # if that category is true then prepend the 0s
                # bc the first 2 cols in the output are not_rare_not_out
                # and not_rare_out
                if out_tb.index == np.array([True]):
                    # case where all gene-ID pairs for the gene
                    # are rare variants
                    out_list = [0] + out_list
                elif len(out_tb.columns) == 1:
                    # case where only gene-ID pair for the gene
                    # is an outlier with a common variant
                    out_list = [0] + out_list
                    out_list.extend([0, 0])
                else:
                    # case where all gene-ID pairs for the gene
                    # are common variants
                    out_list.append(0)
            else:
                # case where all gene-ID pairs are outliers and there are
                # non-zero subsets for both rare and non_rare categories.
                # since the requirement is that we only look at
                # genes with outliers with variants, there must always
                # be a expr_outlier or expr_outlier_neg True column
                out_list.append(0)
                out_list.insert(2, 0)
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
        outlier_df = self.identify_rows_to_keep(
            enrich_df, max_intrapop_af,
            distribution=self.distribution, cut_off_tuple=cut_off_tuple)
        # only keep outliers with rare variants
        outlier_df = outlier_df.loc[
            outlier_df.rare_variant_status & outlier_df.expr_outlier]
        # cols_to_keep = ["blinded_id", "gene", "z_expr", "tss_dist",
        #                 "var_id", "popmax_af", "var_id_freq"]
        # keep only a specific annotation
        # in_anno = outlier_df.loc[:, 'nkx2.5.mm9.hg19'] == 1
        # outlier_df = outlier_df.loc[in_anno]
        # outlier_df = outlier_df[cols_to_keep]
        outlier_df.rename(columns={"var_id_freq": "intra_cohort_af"},
                          inplace=True)
        outlier_df.to_csv(self.rv_outlier_loc, index=False, sep="\t")

#
#
#
#
