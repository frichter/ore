"""Script to count number of annotation bases within 10 kb of every TSS.

module load bedtools/2.27.0
module load python/3.5.0
module load py_packages/3.5
cd /sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_04/
python

"""

import subprocess
import os
import re
import glob
from functools import partial
import multiprocessing as mp

import pandas as pd
# import numpy as np

import pybedtools
from pybedtools import BedTool


pybedtools.set_tempdir('pybedtools_temp_dir/')


"""
# get all genes into a single bed file
cd /sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_04/

cat wgs_atrial_per_chrom/*gene_bed*.bed | wc -l
cat wgs_vent_per_chrom/*gene_bed*.bed | wc -l
cat wgs_arterial_valve_per_chrom/*gene_bed*.bed | wc -l
cat wgs_atrial_per_chrom/*gene_bed*.bed wgs_vent_per_chrom/*gene_bed*.bed \
wgs_arterial_valve_per_chrom/*gene_bed*.bed | wc -l
## 13884+13047+14915 same length
cat wgs_atrial_per_chrom/*gene_bed*.bed wgs_vent_per_chrom/*gene_bed*.bed \
wgs_arterial_valve_per_chrom/*gene_bed*.bed | cut -f1-4 | \
sort -V -k1,1 -k2,2 | uniq > all_genes.bed
"""


"""

cd /hpc/users/richtf01/chdiTrios/Felix/wgs/bed_annotations/hg19_all/
# combine TF binding site annotations
ls *ata4*
# confirm no headers
head -n1 *ata4*
cat *ata4* | sort -V -k1,1 -k2,2 | cut -f1-3 | bedtools merge > any_gata4.bed

ls *kx2*
# confirm no headers
head -n1 *kx2*
cat *kx2* | sort -V -k1,1 -k2,2 | cut -f1-3 | bedtools merge > any_nkx25.bed

ls *bx5*
# confirm no headers
head -n1 *bx5*
cat *bx5* | sort -V -k1,1 -k2,2 | cut -f1-3 | bedtools merge > any_tbx5.bed

cat any_gata4.bed any_nkx25.bed any_tbx5.bed | sort -V -k1,1 -k2,2 | \
    bedtools merge > all_tf.bed

"""


# various file locations
data_dir = ('/hpc/users/richtf01/chdiTrios/Felix/dna_rna/ore/' +
            'ore/data/hg19_genome_masks/')
segdup = data_dir + "hg19_segdup.bed.gz"
lcr = data_dir + "hg19_lcr_hs37d5.bed.gz"
anno_dir = "/hpc/users/richtf01/chdiTrios/Felix/wgs/bed_annotations/hg19_all/"

rmsk = anno_dir + "rmsk.merged.sorted.bed"
map300 = anno_dir + "mappability1_300.bed"
hla_muc = anno_dir + "genes.MUC.HLA.bed"
dac_blacklist = anno_dir + "dac_blacklist.bed"
duke_blacklist = anno_dir + "encode_duke_blacklist.bed"
pseudoauto_XY = anno_dir + "pseudoautosomal_XY.bed"
anno_hg19 = ("/hpc/users/richtf01/chdiTrios/Felix/wgs/bed_annotations/" +
             "hg19_all/*.bed")

gene_bed_f = ("/sc/orga/projects/chdiTrios/Felix/dna_rna/" +
              "wgs_pcgc_2018_04/all_genes.bed")


"""Removing mapping features that were excluded in the outlier analysis.
"""

gene_bed_obj = BedTool(gene_bed_f)

print(gene_bed_obj.head())
print(gene_bed_obj.slop(genome="hg19", b=1e4).head())
print(gene_bed_obj.slop(genome="hg19", b=1e4).count())
print(gene_bed_obj.count())

# add 10 kb to flanks
gene_bed_obj = gene_bed_obj.slop(genome="hg19", b=1e4)

file_loc_list = [rmsk, segdup, lcr, map300, hla_muc, dac_blacklist,
                 duke_blacklist, pseudoauto_XY]

rm_list = []
count = 0
for file_loc in file_loc_list:
    bed_iterable = glob.iglob(file_loc)
    for bed_name in bed_iterable:
        print(bed_name)
        rm_list.append(bed_name)
        bed = BedTool(bed_name)
        gene_bed_obj = gene_bed_obj.subtract(bed)
        count += 1
        print(count)


print(gene_bed_obj.head())
print(gene_bed_obj.count())
gene_bed_obj.saveas(re.sub(".bed", "_hi_mapping.bed", gene_bed_f))


"""Overlapping the gene TSS bed file with each sorted/merged annotation.
"""


def overlap_gene_bed(bed_name, gene_bed_obj, out_dir):
    """Overlap annotations with 10 kbflanked gene TSSs."""
    print(bed_name)
    if bed_name.endswith("human_acc_regions_heart_enh.bed"):
        return "Never run: " + bed_name
    out_f = re.sub(".*/", out_dir, bed_name)
    sorted_f = re.sub(".bed$", "_sorted.bed", out_f)
    ct_out_f = re.sub(".bed$", "_gene_ct.bed", out_f)
    if os.path.exists(ct_out_f):
        return "Not rerun: " + bed_name
    # sort the bed file (necessary for merging)
    subprocess.call('sort -V -k1,1 -k2,2 {0} > {1}'.format(bed_name, sorted_f),
                    shell=True)
    bed = BedTool(sorted_f)
    # merge the bed file before counting the intersection
    bed = bed.merge()
    gene_bed_obj = gene_bed_obj.intersect(bed, wao=True)
    gene_bed_obj.saveas(ct_out_f)
    return bed_name


# cd /hpc/users/richtf01/chdiTrios/Felix/wgs/bed_annotations/\
# gene_feature_count_hg19/
out_dir = ("/hpc/users/richtf01/chdiTrios/Felix/wgs/bed_annotations/" +
           "gene_feature_count_hg19/")
gene_bed_f = ("/sc/orga/projects/chdiTrios/Felix/dna_rna/" +
              "wgs_pcgc_2018_04/all_genes.bed")
gene_bed_obj = BedTool(re.sub(".bed", "_hi_mapping.bed", gene_bed_f))
bed_iterable = glob.iglob(anno_hg19)
partial_overlap_gene_bed = partial(
    overlap_gene_bed, gene_bed_obj=gene_bed_obj, out_dir=out_dir)
# testing or running sequentially
count = 0
anno_completed = []
for bed_name in bed_iterable:
    bed_done = partial_overlap_gene_bed(bed_name)
    anno_completed.append(bed_done)
    count += 1
    # if count > 3:
    #     break

[i for i in anno_completed if not i.startswith("Not")]
pool = mp.Pool(processes=6)
# print("Total available cores: " + str(mp.cpu_count()))
anno_completed = pool.map(partial_overlap_gene_bed, bed_iterable)

len(anno_completed)

"""Load and clean all overlap counts into a single dataframe/file.
"""


def load_and_clean_cts(bed_file):
    """Load BED overlap data and clean column names."""
    # clean name
    rep_w_blank = (".*/|.merged.sorted|.sorted|.bed$|.bed.gz$|.txt$|_gene_ct" +
                   "|_sorted")
    bed_name = re.sub(rep_w_blank, "", bed_file)
    bed_name = re.sub("all_predictions", "cvdc_enhancers_dickel", bed_name)
    # use as column name
    print(bed_name)
    bed_df = pd.read_table(bed_file, header=None, index_col=False)
    # names=bed_header
    bed_col_count = bed_df.shape[1]
    bed_header = ["Chrom", "Start", "End", "Gene", "overlap_ct"]
    extra_coL_ct = 0
    # get rid of all the extra/useless columns
    while len(bed_header) < bed_col_count:
        extra_coL_ct += 1
        bed_header.insert(4, "rm" + str(extra_coL_ct))
    bed_df.columns = bed_header
    if not bed_df.Chrom.values[0].startswith("chr"):
        print(bed_df.head())
        raise ImportError()
    bed_df = bed_df[["Gene", "overlap_ct"]]
    bed_df["anno"] = bed_name
    bed_df = bed_df.groupby(["Gene", "anno"]).sum()
    return bed_df


out_dir = ("/hpc/users/richtf01/chdiTrios/Felix/wgs/bed_annotations/" +
           "gene_feature_count_hg19/")
bed_iterable = glob.iglob(out_dir + "*_gene_ct.bed")

# For Testing:
# count = 0
# for bed_file_i in bed_iterable:
#     bed_df_i = load_and_clean_cts(bed_file_i)
#     print(bed_df_i.head())
#     print(bed_df_i.shape)
#     count += 1
#     if count > 10:
#         break


partial_load_and_clean_cts = partial(load_and_clean_cts)
# keep partial code in case you want to add other arguments in future
pool = mp.Pool(processes=6)
# print("Total available cores: " + str(mp.cpu_count()))
bed_df_list = pool.map(partial_load_and_clean_cts, bed_iterable)
[bed_i.shape for bed_i in bed_df_list]
len(bed_df_list)
bed_df = pd.concat(bed_df_list)

bed_df.columns
bed_df.head()
bed_df.shape


gene_bed_f = ("/sc/orga/projects/chdiTrios/Felix/dna_rna/" +
              "wgs_pcgc_2018_04/all_genes.bed")
# saving
bed_df.to_csv(re.sub(".bed", "_anno_cts_grouped_from_merged.txt", gene_bed_f),
              sep="\t")  # keep index after groupby so don't use index=False
# _anno_cts.txt _anno_cts_grouped.txt
# mp result shape: (26624244, 1)
# non-mp results shape: (26624242, 1)
