#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Create dictionary with WES ID (keys): WGS IDs (values) from QC data

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-07-07
:Copyright: 2018, Felix Richter
:License: CC BY-SA

cd /Users/felixrichter/Dropbox/PhD/alzheimers/metadata
python3
"""

import pandas as pd

mt_dir = '/Users/felixrichter/Dropbox/PhD/alzheimers/metadata/'

# load WGS IDs
with open(mt_dir + 'wgs_samples_2018_05.txt', 'r') as f:
    wgs_id_list = ['wgs' + i.strip() for i in f]


# load QC dataframe
cols_to_keep = ['Sampleid', 'Datatype', 'BrainID by annotation',
                'BrainID (final)', 'Action']
qc_df = pd.read_table(mt_dir + 'MSBB_RNAseq.WES.WGS_sample_QC_info.csv',
                      sep=',', usecols=cols_to_keep)
qc_df.head()
qc_df.shape
len(qc_df.Sampleid.unique())
# all 349 samples are represented (I hope)
len([i for i in wgs_id_list if any(qc_df.Sampleid.str.contains(i))])

wgs_line = qc_df.Sampleid.isin(wgs_id_list)
keep_line = qc_df.Action.str.contains('OKay')
wgs_df = qc_df.loc[wgs_line & keep_line]
de_id_list = wgs_df['BrainID (final)'].tolist()

wes_df = qc_df.loc[qc_df.Datatype == 'WES']
mismatch_wes_id_df = wes_df[
    wes_df['BrainID by annotation'] != wes_df['BrainID (final)']]

# do the mismatched IDs have RNAseq?
rna_df = qc_df.loc[qc_df.Datatype == 'RNA-seq']

mismatch_w_rnaseq = rna_df['BrainID (final)'].isin(
    mismatch_wes_id_df['BrainID by annotation'])
rna_df[mismatch_w_rnaseq]
mismatch_w_rnaseq = rna_df['BrainID (final)'].isin(
    mismatch_wes_id_df['BrainID (final)'])
rna_df[mismatch_w_rnaseq]


final_col = wes_df['BrainID (final)'].isin(de_id_list)
initial_col = wes_df['BrainID by annotation'].isin(de_id_list)
wes_okay = wes_df.Action.str.contains('OKay')
wes_df = wes_df[final_col & initial_col & wes_okay]

wgs_df = wgs_df.loc[:, ['Sampleid', 'BrainID (final)']]
wgs_df.set_index('BrainID (final)', inplace=True)
wgs_df.rename(columns={"Sampleid": "wgs_id"}, inplace=True)
wes_df = wes_df.loc[:, ['Sampleid', 'BrainID (final)']]
wes_df.set_index('BrainID (final)', inplace=True)
wes_df.rename(columns={"Sampleid": "wes_id"}, inplace=True)
final_mapping_df = wgs_df.join(wes_df, how='inner')
final_mapping_df.wgs_id = final_mapping_df.wgs_id.str.replace('wgs', '')
final_mapping_df.head()
final_mapping_df.shape

wgs_wes_id_loc = mt_dir + 'id_map_wgs_wes_2018_08_22.txt'
final_mapping_df.to_csv(wgs_wes_id_loc, sep='\t', index=False)

# wes ID should be key, wgs ID should be value
wgs_wes_dict = {}

#
