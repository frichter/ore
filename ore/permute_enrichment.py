#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Calculate permutation p-values for enrichment.

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-09-13
:Copyright: 2018, Felix Richter
:License: CC BY-SA
"""

import copy
import glob
import re
import gc
from datetime import datetime
import csv

import pandas as pd
import numpy as np

from .enrichment import Enrich
from .utils import prepare_directory


class PermuteEnrich(Enrich):
    """Permute inputs to enrichment functions."""

    def __init__(self, joined_df, expr_outlier_df, output_prefix,
                 distribution, anno_list, obs_enrich_loc,
                 loop_enrich_args, write_rv_args, n_perms=1):
        """Initialize a permuted enrichment object."""
        self.joined_df = joined_df
        self.drop_relevant_expression_columns()
        self.expr_outlier_df = expr_outlier_df
        # other inputs
        self.distribution = distribution
        self.anno_list = anno_list
        self.obs_enrich_loc = re.sub(".txt$", "_gene.txt", obs_enrich_loc)
        perm_dir = output_prefix + '_per_chrom/perms'
        prepare_directory(perm_dir)  # , clean_run=True
        self.generic_enrich_loc = perm_dir + '/perm_{}_{}.txt'
        self.loop_enrich_args = loop_enrich_args
        self.write_rv_args = write_rv_args
        # loop over permutations
        for n_perm in range(n_perms):
            print('=========== Permutation {} ==========='.format(str(n_perm)))
            self.run_permutation(n_perm)
            gc.collect()
        self.compare_observed_to_permuted()

    def drop_relevant_expression_columns(self):
        """Remove expression-relevant columns that will be permuted."""
        cols_to_drop = ['z_expr', 'expr_rank', 'z_abs', 'expr_outlier_neg',
                        'expr_outlier_pos', 'expr_outlier']
        cols_to_drop = [i for i in cols_to_drop if i in self.joined_df.columns]
        self.joined_df.drop(cols_to_drop, axis=1, inplace=True)

    def run_permutation(self, n_perm):
        """Run a single permutation."""
        # randomly re-assign IDs in outlier_df
        self.permute_ids()
        # join outlier_df_permute with all_data
        self.joined_df.set_index(['gene', 'blinded_id'], inplace=True)
        self.joined_df = self.joined_df.join(self.permute_expr_df, how='inner')
        # declare enrichment output file name
        ts = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
        self.enrich_loc = self.generic_enrich_loc.format(str(n_perm), ts)
        self.rv_outlier_loc = re.sub('.txt$', '_rv_outs.txt', self.enrich_loc)
        # run loop_enrichment
        self.joined_df.reset_index(inplace=True)
        print(self.joined_df.head())
        super(PermuteEnrich, self).write_rvs_w_outs_to_file(
            **self.write_rv_args)
        super(PermuteEnrich, self).loop_enrichment(**self.loop_enrich_args)
        # (possibly just a faster method to identify the number of
        #  outliers with rare variants)
        self.drop_relevant_expression_columns()

    def permute_ids(self):
        """Permute the RNAseq IDs."""
        uniq_ids = self.expr_outlier_df.index.get_level_values(
            'blinded_id').unique().tolist()
        # using permutation instead of shuffle:
        # https://stackoverflow.com/a/15474335
        perm_ids = np.random.permutation(uniq_ids).tolist()
        self.id_dict = dict(zip(uniq_ids, perm_ids))
        print(self.id_dict)
        self.permute_expr_df = copy.deepcopy(self.expr_outlier_df)
        self.permute_expr_df.reset_index(inplace=True)
        # print(self.permute_expr_df.head())
        self.permute_expr_df = self.permute_expr_df.assign(
            blinded_id=self.permute_expr_df['blinded_id'].map(self.id_dict))
        self.permute_expr_df.set_index(['gene', 'blinded_id'], inplace=True)
        # print(self.permute_expr_df.head())

    def compare_observed_to_permuted(self):
        """Compare results from observed and permuted data to get p-values."""
        obs_df = pd.read_table(self.obs_enrich_loc)
        print("Observed:")
        (n_nom_sig_obs, max_or_obs, min_p_obs,
            rv_outs_obs) = self.get_sig_metrics(obs_df)
        print(n_nom_sig_obs, max_or_obs, min_p_obs, rv_outs_obs)
        # load permutation enrichments 1 at a time
        perm_f_iter = glob.iglob(self.generic_enrich_loc.format('*', '*_gene'))
        (total_perms, total_nom_sig, total_max_or, total_min_p,
            total_rv_outs) = 0, 0, 0, 0, 0
        perm_dict = {'perm_ct': [], 'nom_sig_list': [],
                     'max_or_list': [], 'min_p_list': [],
                     'rv_outs_list': []}
        print("Reviewing permutations:")
        for perm_f in perm_f_iter:
            # print(perm_f)
            perm_df = pd.read_table(perm_f)
            n_nom_sig, max_or, min_p, rv_outs = self.get_sig_metrics(perm_df)
            total_nom_sig += int(n_nom_sig >= n_nom_sig_obs)
            total_max_or += int(max_or >= max_or_obs)
            total_min_p += int(min_p <= min_p_obs)
            total_rv_outs += int(rv_outs >= rv_outs_obs)
            total_perms += 1
            perm_dict['perm_ct'].append(total_perms)
            perm_dict['nom_sig_list'].append(n_nom_sig)
            perm_dict['max_or_list'].append(max_or)
            perm_dict['min_p_list'].append(min_p)
            perm_dict['rv_outs_list'].append(rv_outs)
            if total_perms % 10 == 0:
                print(total_perms, total_nom_sig, total_max_or, total_min_p,
                      total_rv_outs)
            # if total_perms > 999:
            #     break
        print("Permutation p-values are:")
        print("Nominally sig: " + str(total_nom_sig/total_perms))
        print("Max OR: " + str(total_max_or/total_perms))
        print("Min p: " + str(total_min_p/total_perms))
        print("RV out count: " + str(total_rv_outs/total_perms))
        print(n_nom_sig_obs, max_or_obs, min_p_obs, rv_outs_obs)
        ts = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
        perm_stats_f = self.generic_enrich_loc.format(
            'summary', 'stats_' + ts)
        with open(perm_stats_f, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(perm_dict.keys())
            writer.writerows(zip(*perm_dict.values()))
        # run check_permutation_significance
        # count number of examples as or more extreme
        # print a histogram of the results

    def get_sig_metrics(self, df):
        """Test if a permutation is more extreme for multiple metrics."""
        # identify if permutation is more extreme based on:
        # if total number of RV-outlier pairs more than observed. (Gabriel)
        # if number of nominally significant associations (P<0.05, OR>1)
        df = df[df.tss_cut_off == 1e4]
        rv_outs = df[df.af_cut_off == 1e-5]['rare_out'].values[0]
        # df_nom_sig = df[(df.p < 0.05) & (df['or'] > 1)]
        # just using those with OR>1 for 1-sided tests
        df_nom_sig = df[df['or'] > 1]
        # 2-sided use:
        # df_nom_sig = df
        n_nom_sig = df_nom_sig[df_nom_sig.p < 0.05].shape[0]
        if df_nom_sig.shape[0] == 0:
            max_or = 0
            min_p = 1
        else:
            # if most extreme permutation is more than most extreme observed
            max_or = df_nom_sig['or'].max()
            min_p = df_nom_sig.p.min()
            # number of RV outs is more for most significant observed
            # max_or_line = df_nom_sig['or'] == max_or
        return n_nom_sig, max_or, min_p, rv_outs


#
#
#


#
