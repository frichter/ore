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
            # ask if header ID among keys of gt_to_wgs_id_dict
            print(len([i for i in header if i in gt_to_wgs_id_dict]))
            line_indices_to_keep.extend([i for i, gt_id in enumerate(header)
                                         if gt_id in gt_to_wgs_id_dict])
            # rna_id_to_keep = [rna_id for rna_id in header
            #                   if rna_id in rna_data.id_dict]
            new_header = header[0:4]
            new_header.extend([gt_to_wgs_id_dict[gt_id] for gt_id in header
                               if gt_id in gt_to_wgs_id_dict])
            # print(new_header)
            print("Confirm NEW header columns are unique " +
                  "(that these numbers are the same):")
            print(len(new_header), len(set(new_header)))
            with open(rna_data.loc_w_wgs_id, 'w') as expr_w_wgs_id:
                _ = expr_w_wgs_id.write("\t".join(new_header) + "\n")
                for line in expr_f:
                    line_list = str(line, 'utf-8').split("\t")
                    line_to_keep = [line_i for i, line_i in
                                    enumerate(line_list) if i in
                                    line_indices_to_keep]
                    out_line = "\t".join(line_to_keep)  # + "\n"
                    _ = expr_w_wgs_id.write(out_line)
                    _


"""Load dictionary mapping WGS IDs (values) to WES ID (keys)."""
mt_loc = ("/hpc/users/richtf01/alzheimers/metadata/" +
          "id_map_wgs_wes_2018_08_22.txt")
gt_to_wgs_id_dict = {}
with open(mt_loc) as f:
    header = next(f)
    for line in f:
        val, key = line.split()
        gt_to_wgs_id_dict[key] = val

len(gt_to_wgs_id_dict)

"""Updated the header per RNAseq experiment."""
# expr_loc = ("residuals_AMPAD_MSSM_GE_SV_17_tissue_44_with_disease_in_" +
#             "model_europeans_only.bed.gz") _new_z
expr_f_list = [i for i in glob.iglob("residuals_AMPAD*_only.bed.gz")]

# expr_f_list[0] expr_f_list[1]
for expr_f in expr_f_list:
    rna_data = RNAdata(expr_f)
    rna_data.update_header(gt_to_wgs_id_dict)


"""bgzip and tabix RNAseq experiments with WGS ID headers."""
for expr_f in expr_f_list:
    rna_data = RNAdata(expr_f)
    # bgzip and re-index new file
    tbx_cmd = "time bgzip {} && time tabix -p bed {}.gz".format(
        rna_data.loc_w_wgs_id, rna_data.loc_w_wgs_id)
    print(tbx_cmd)
    subprocess.call(tbx_cmd, shell=True)


"""
NEW and improved final kept:
10: 175/188
22: 158/172
44: 148/158
36: 135/147

"""


#
#
#
