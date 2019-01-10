#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Utilities for calculating enrichment.

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-05-20
:Copyright: 2018, Felix Richter
:License: CC BY-SA
"""


import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def calculate_var_enrichment(enrich_df, expr_df):
    """Calculate variant-centric enrichment."""
    out_genes = expr_df[expr_df.expr_outlier].index.get_level_values(
        'gene').unique()
    enrich_df = enrich_df[enrich_df.gene.isin(out_genes)]
    out_tb = pd.crosstab(enrich_df.rare_variant_status,
                         enrich_df.expr_outlier)
    # print(out_tb)
    out_list = flatten_crosstab(out_tb)
    # only keep subset of genes that are negative expression outliers
    out_genes = expr_df[expr_df.expr_outlier_neg].index.get_level_values(
        'gene').unique()
    enrich_df_neg = enrich_df[enrich_df.gene.isin(out_genes)]
    neg_out_tb = pd.crosstab(enrich_df_neg.rare_variant_status,
                             enrich_df_neg.expr_outlier_neg)
    # print(neg_out_tb)
    neg_out_list = flatten_crosstab(neg_out_tb)
    # now positive outliers
    out_genes = expr_df[expr_df.expr_outlier_pos].index.get_level_values(
        'gene').unique()
    enrich_df_pos = enrich_df[enrich_df.gene.isin(out_genes)]
    pos_out_tb = pd.crosstab(enrich_df_pos.rare_variant_status,
                             enrich_df_pos.expr_outlier_pos)
    # print(pos_out_tb)
    pos_out_list = flatten_crosstab(pos_out_tb)
    # now perform the actual calculations (if possible)
    try:
        fet_or, fet_p = fisher_exact(out_tb)
        ci_lo, ci_hi = calculate_ci(fet_or, out_list)
    except ValueError:
        fet_or, fet_p, ci_lo, ci_hi = ("NA", "NA", "NA", "NA")
    try:
        fet_or_neg, fet_p_neg = fisher_exact(neg_out_tb)
        ci_neg_lo, ci_neg_hi = calculate_ci(fet_or_neg, neg_out_list)
    except ValueError:
        fet_or_neg, fet_p_neg, ci_neg_lo, ci_neg_hi = ("NA", "NA", "NA", "NA")
    try:
        fet_or_pos, fet_p_pos = fisher_exact(pos_out_tb)
        ci_pos_lo, ci_pos_hi = calculate_ci(fet_or_pos, pos_out_list)
    except ValueError:
        fet_or_pos, fet_p_pos, ci_pos_lo, ci_pos_hi = ("NA", "NA", "NA", "NA")
    out_list.extend(neg_out_list)
    out_list.extend(pos_out_list)
    out_list.extend([fet_or, fet_p, ci_lo, ci_hi, fet_or_neg,
                     fet_p_neg, ci_neg_lo, ci_neg_hi,
                     fet_or_pos, fet_p_pos, ci_pos_lo, ci_pos_hi])
    print(out_list)
    return out_list


def calculate_gene_enrichment(enrich_df, expr_df):
    """Calculate gene-centric enrichment.

    Counts each gene-ID pair once (i.e., removes extra variants per
        gene-ID pair with `drop_duplicates`)

    """
    enrich_df['gene_has_rare_var'] = enrich_df.groupby(
        ['gene', 'blinded_id'])['rare_variant_status'].transform('sum') > 0
    enrich_df = enrich_df.loc[:, [
        'gene', 'blinded_id', 'expr_outlier', 'expr_outlier_neg',
        'gene_has_NEG_out_w_vars', 'expr_outlier_pos',
        'gene_has_POS_out_w_vars',
        'gene_has_rare_var']
        ].drop_duplicates(keep='first')
    # count number of outliers and non-outliers with gene_has_rare_var
    out_class = 'expr_outlier'
    out_list, out_tb = quantify_gene_outs(out_class, enrich_df, expr_df)
    # only keep subset of genes that are negative expression outliers
    out_class = 'expr_outlier_neg'
    neg_out_list, neg_out_tb = quantify_gene_outs(
        out_class, enrich_df, expr_df)
    # now positive outliers
    out_class = 'expr_outlier_pos'
    pos_out_list, pos_out_tb = quantify_gene_outs(
        out_class, enrich_df, expr_df)
    # now perform the actual calculations (if possible)
    try:
        fet_or, fet_p = fisher_exact(out_tb)
        ci_lo, ci_hi = calculate_ci(fet_or, out_list)
    except ValueError:
        fet_or, fet_p, ci_lo, ci_hi = ("NA", "NA", "NA", "NA")
    try:
        fet_or_neg, fet_p_neg = fisher_exact(neg_out_tb)
        ci_neg_lo, ci_neg_hi = calculate_ci(fet_or_neg, neg_out_list)
    except ValueError:
        fet_or_neg, fet_p_neg, ci_neg_lo, ci_neg_hi = ("NA", "NA", "NA", "NA")
    try:
        fet_or_pos, fet_p_pos = fisher_exact(pos_out_tb)
        ci_pos_lo, ci_pos_hi = calculate_ci(fet_or_pos, pos_out_list)
    except ValueError:
        fet_or_pos, fet_p_pos, ci_pos_lo, ci_pos_hi = ("NA", "NA", "NA", "NA")
    out_list.extend(neg_out_list)
    out_list.extend(pos_out_list)
    out_list.extend([fet_or, fet_p, ci_lo, ci_hi, fet_or_neg,
                     fet_p_neg, ci_neg_lo, ci_neg_hi,
                     fet_or_pos, fet_p_pos, ci_pos_lo, ci_pos_hi])
    print(out_list)
    return out_list


def quantify_gene_outs(out_class, enrich_df, expr_df):
    """Create the 2x2 table for expression outliers."""
    out_genes = expr_df[expr_df[out_class]].index.get_level_values(
        'gene').unique()
    enrich_df = enrich_df[enrich_df.gene.isin(out_genes)]
    rare_out = enrich_df[enrich_df.gene_has_rare_var & enrich_df[out_class]
                         ].shape[0]
    # print("All genes:")
    # print(enrich_df.shape)
    enrich_df_w_outs = enrich_df[enrich_df.gene.isin(out_genes)]
    # print("Only genes with " + out_class + ":")
    # print(enrich_df_w_outs.shape)
    # print("Unique outlier genes: " + str(len(out_genes)))
    rare_not_out = enrich_df_w_outs[enrich_df_w_outs.gene_has_rare_var &
                                    (~enrich_df_w_outs[out_class])].shape[0]
    not_rare_out = expr_df[out_class].sum() - rare_out
    # count number of outliers in long expression DF
    # print(out_class + " unique outlier genes: " + str(len(out_genes)))
    total_gene_by_id = expr_df.loc[expr_df.index.get_level_values(
        'gene').isin(out_genes)].shape[0]
    # subtract sum of list from the number of rows in long expression DF
    not_rare_not_out = (
        total_gene_by_id - (rare_out + rare_not_out + not_rare_out))
    out_list = [not_rare_not_out, not_rare_out, rare_not_out, rare_out]
    out_tb = np.array([[not_rare_not_out, not_rare_out],
                       [rare_not_out, rare_out]])
    return out_list, out_tb


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


def calculate_ci(odds_ratio, val_list):
    """Calculate confidence intervals from FET results.

    Sources
        https://stats.stackexchange.com/a/1483
        https://stats.stackexchange.com/a/2233

    Warnings to address
        RuntimeWarning: invalid value encountered in double_scalars
        RuntimeWarning: divide by zero encountered in log
    """
    val_array = np.array(val_list).astype(float)
    or_se = sum(np.reciprocal(val_array))**(1/2.0)
    ci = (np.exp(np.log(odds_ratio) - 1.96*or_se),
          np.exp(np.log(odds_ratio) + 1.96*or_se))
    return ci


#
#
#
