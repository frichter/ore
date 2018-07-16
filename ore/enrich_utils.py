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
    out_list = flatten_crosstab(out_tb)
    # out_list = out_tb.values.flatten().tolist()
    # while len(out_list) < 4:
    #     out_list.append("NA")
    # only keep subset of genes with low expression outliers
    enrich_df = enrich_df[enrich_df.gene_has_NEG_out_w_vars]
    neg_out_tb = pd.crosstab(enrich_df.gene_has_rare_var,
                             enrich_df.expr_outlier_neg)
    print(neg_out_tb)
    neg_out_list = flatten_crosstab(neg_out_tb)
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
