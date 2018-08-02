#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Update headers with correct ID

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-07-07
:Copyright: 2018, Felix Richter
:License: CC BY-SA

cd /sc/orga/projects/chdiTrios/Felix/alzheimers/expression
module load samtools/1.3
module load python/3.5.0
module load py_packages/3.5
python

"""

import re
import subprocess
import gzip
import glob


def load_id_map_dict(tissue, id_dict):
    """Load map between WGS and Sample IDs."""
    id_map_loc = ("/hpc/users/richtf01/alzheimers/metadata/wgs_idt" +
                  "ornasq_id/msbb.BM_{}.PMI_AOD_race_sex_RIN_" +
                  "exonicRate_rRnaRate_batch_adj.WGSoverlapped" +
                  ".SampleOrder").format(tissue)
    with open(id_map_loc, 'r') as id_map:
        header = next(id_map)
        print(header)
        for line in id_map:
            line_list = line.strip().split("\t")
            rna_id = re.sub(".*BM", "BM", line_list[2])
            rna_id = re.sub(".*hB_RNA", "hB_RNA", rna_id)
            sample_id = line_list[1]
            if (sample_id in id_dict) and (id_dict[sample_id] !=
                                           line_list[0]):
                print("{} already in dict for {}, confirm same as {}".format(
                    sample_id, id_dict[sample_id], line_list[0]))
            id_dict[sample_id] = line_list[0]
    return id_dict


sample_to_wgs_id_dict = {}
for tissue in ["10", "22", "36", "44"]:
    print(tissue)
    sample_to_wgs_id_dict = load_id_map_dict(tissue, sample_to_wgs_id_dict)


# got 243 samples with wgs_ids. However I have 350 WGS ids..
len(sample_to_wgs_id_dict)

gt_to_wgs_id_dict = {}
mt_loc = "/hpc/users/richtf01/alzheimers/metadata/meta.BM36.txt"
with open(mt_loc, 'r') as mt:
    header = next(mt)
    for line in mt:
        line_list = line.strip().split("\t")
        sample_id = line_list[0]
        gt_id = re.sub("-", "_", line_list[2]) + "_" + line_list[1]
        print(sample_id)
        print(gt_id)
        if gt_id in gt_to_wgs_id_dict:
            print("{} already in dict for {}, confirm same as {}".format(
                gt_id, gt_to_wgs_id_dict[gt_id],
                sample_to_wgs_id_dict[sample_id]))
        if sample_id in sample_to_wgs_id_dict:
            gt_to_wgs_id_dict[gt_id] = sample_to_wgs_id_dict[sample_id]


# only have 161 genotype IDs...
len(gt_to_wgs_id_dict)


class RNAdata(object):
    """Methods and objects for manipulating RNAseq."""
    
    def __init__(self, expr_loc):
        """Initialize RNAseq object."""
        self.gz_loc = expr_loc
        self.loc = expr_loc[:-3]
        self.loc_w_wgs_id = re.sub(".bed.gz", "_wgs_ids.bed", expr_loc)
        self.tissue = re.sub(".*_tissue_", "", self.loc)
        self.tissue = re.sub("_with_disease_in.*", "", self.tissue)
        # self.id_map_loc = ("/hpc/users/richtf01/alzheimers/metadata/wgs_idt"+
        #                    "ornasq_id/msbb.BM_{}.PMI_AOD_race_sex_RIN_" +
        #                    "exonicRate_rRnaRate_batch_adj.WGSoverlapped" +
        #                    ".SampleOrder").format(self.tissue)
        # self.sample_to_wgs_id_dict = sample_to_wgs_id_dict
        # self.gt_to_sample_id_dict = gt_to_sample_id_dict
        print(self.loc)
        print(self.loc_w_wgs_id)
        print(self.tissue)
    
    def update_header(self, gt_to_wgs_id_dict):
        """Update the header and keep only data for samples with a WGS id."""
        with gzip.open(rna_data.gz_loc, 'r') as expr_f:
            header = str(next(expr_f), 'utf-8').strip().split("\t")
            # which rnaseq IDs have a WGS id?
            header_gt_ids = header[4:]
            line_indices_to_keep = [i for i in range(0, 4)]
            # print(header_rna_ids)
            # print([i for i in header if i in gt_to_wgs_id_dict])
            print("Total IDs:", len(header_gt_ids))
            print("IDs being kept:")
            print(len([i for i in header if i in gt_to_wgs_id_dict]))
            line_indices_to_keep.extend([i for i, gt_id in enumerate(header)
                                         if gt_id in gt_to_wgs_id_dict])
            # rna_id_to_keep = [rna_id for rna_id in header
            #                   if rna_id in rna_data.id_dict]
            new_header = header[0:4]
            new_header.extend([gt_to_wgs_id_dict[gt_id] for gt_id in header
                               if gt_id in gt_to_wgs_id_dict])
            # print(new_header)
            print("Confirm all IDs unique (that these numbers are the same:)")
            print(len(new_header), len(set(new_header)))
            with open(rna_data.loc_w_wgs_id, 'w') as expr_w_wgs_id:
                _ = expr_w_wgs_id.write("\t".join(new_header) + "\n")
                for line in expr_f:
                    line_list = str(line, 'utf-8').split("\t")
                    line_to_keep = [line_i for i, line_i in
                                    enumerate(line_list) if i in
                                    line_indices_to_keep]
                    out_line = "\t".join(line_to_keep) + "\n"
                    _ = expr_w_wgs_id.write(out_line)
                    _


# expr_loc = ("residuals_AMPAD_MSSM_GE_SV_17_tissue_36_with_disease_in_" +
#             "model_europeans_only.bed.gz")
# expr_loc = ("residuals_AMPAD_MSSM_GE_SV_17_tissue_44_with_disease_in_" +
#             "model_europeans_only.bed.gz")
expr_f_list = [i for i in glob.iglob("residuals_AMPAD*_only_new_z.bed.gz")]

# expr_f_list[0] expr_f_list[1]
for expr_f in expr_f_list:
    rna_data = RNAdata(expr_f)
    rna_data.update_header(gt_to_wgs_id_dict)


# bgzip and re-index new file
tbx_cmd = "time bgzip {} && time tabix -p bed {}.gz".format(
    rna_data.loc_w_wgs_id, rna_data.loc_w_wgs_id)
print(tbx_cmd)
subprocess.call(tbx_cmd, shell=True)

"""
Final kept:
10: 76/188
22: 69/172
44: 69/158
36: 86/147

Sources of ID loss:
Only 243 WGS IDs with BB sample IDs received from Chen
243/348
Only 161/243 of these BB sample IDs have genotype IDs in meta.BM36.txt
len(sample_to_wgs_id_dict)
len(gt_to_wgs_id_dict)

"""

"""
Recalculate z-scores (possibly add as code to ORE)

import glob
import re
import subprocess

import pandas as pd

expr_f_list = [i for i in glob.iglob("residuals_AMPAD*_only.bed.gz")]

expr_f_list_i = expr_f_list[2]
expr_df = pd.read_table(expr_f_list_i, sep='\t')
expr_df.shape
expr_df.head()
expr_df.rename(columns={expr_df.columns[3]: "gene"}, inplace=True)
expr_df.set_index(['gene'], inplace=True)

gene_df = expr_df.iloc[:, :3]
gene_df.head()
expr_df = expr_df.iloc[:, 3:]

expr_df.shape
expr_df_mean = expr_df.mean(axis=1)
expr_df_std = expr_df.std(axis=1)
expr_df = expr_df.sub(expr_df_mean, axis=0)
expr_df = expr_df.div(expr_df_std, axis=0)

expr_df.head()
expr_df.shape

expr_df = gene_df.join(expr_df, how='inner')
expr_df.head()
expr_df.shape

# reorder columns
expr_df.reset_index(inplace=True)
cols = expr_df.columns.tolist()
cols = cols[1:4] + cols[0:1] + cols[4:]
expr_df = expr_df.loc[:, cols]


expr_df_loc = re.sub(".bed.gz", "_new_z.bed", expr_f_list_i)
expr_df.to_csv(expr_df_loc, index=False, sep="\t", float_format='%g')

tbx_cmd = "time bgzip {} && time tabix -p bed {}.gz".format(
    expr_df_loc, expr_df_loc)
print(tbx_cmd)
subprocess.call(tbx_cmd, shell=True)


"""


#
#
#
