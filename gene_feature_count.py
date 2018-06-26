"""Script to count number of annotation bases within 10 kb of every TSS.

module load bedtools/2.27.0
module load python/3.5.0
module load py_packages/3.5
cd /sc/orga/projects/chdiTrios/Felix/dna_rna/wgs_pcgc_2018_04/
python

"""

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


def overlap_gene_bed(bed_name, gene_bed_obj, out_dir):
    """Overlap annotations with 10 kbflanked gene TSSs."""
    print(bed_name)
    out_f = re.sub(".*/", out_dir, bed_name)
    out_f = re.sub(".bed$", "_gene_ct.bed", out_f)
    bed = BedTool(bed_name)
    gene_bed_obj = gene_bed_obj.intersect(bed, wao=True)
    gene_bed_obj.saveas(out_f)
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
pool = mp.Pool(processes=6)
# print("Total available cores: " + str(mp.cpu_count()))
anno_completed = pool.map(partial_overlap_gene_bed, bed_iterable)


"""
file_loc_list = [anno_hg19]

overlap_list = []
count = 0
for file_loc in file_loc_list:
    bed_iterable = glob.iglob(file_loc)
    for bed_name in bed_iterable:
        print(bed_name)
        out_f = re.sub(".*/", out_dir, bed_name)
        out_f = re.sub(".bed$", "_gene_ct.bed", bed_name)
        overlap_list.append(bed_name)
        bed = BedTool(bed_name)
        gene_bed_obj = gene_bed_obj.intersect(bed, wao=True)
        count += 1
        print(count)
        # if count > 5:
        #     break

print(gene_bed_obj.head())
print(gene_bed_obj.count())
"""

"""Load and clean all overlap counts into a single dataframe/file.
"""


def load_and_clean_cts(bed_file, bed_header):
    """Load BED overlap data and clean column names."""
    # clean name
    rep_w_blank = ".*/|.merged.sorted|.sorted|.bed$|.bed.gz$|.txt$|_gene_ct"
    bed_name = re.sub(rep_w_blank, "", bed_file)
    bed_name = re.sub("all_predictions", "cvdc_enhancers_dickel", bed_name)
    # use as column name
    print(bed_name)
    bed_df = pd.read_table(bed_file, names=bed_header)
    bed_df = bed_df[["Gene", "overlap_ct"]]
    bed_df["anno"] = bed_name
    bed_df = bed_df.groupby(["Gene", "anno"]).sum()
    return bed_df


bed_header = ["Chrom", "Start", "End", "Gene", "rm1", "rm2", "rm3",
              "overlap_ct"]
bed_iterable = glob.iglob(out_dir + "*_gene_ct.bed")

partial_load_and_clean_cts = partial(load_and_clean_cts, bed_header=bed_header)
pool = mp.Pool(processes=6)
# print("Total available cores: " + str(mp.cpu_count()))
bed_df_list = pool.map(partial_load_and_clean_cts, bed_iterable)
len(bed_df_list)
bed_df = pd.concat(bed_df_list)

# bed_df = partial_load_and_clean_cts([i for i in bed_iterable][0])


cols_to_keep = ["Gene", "overlap_ct"]
# len([i for i in bed_iterable])
bed_df_list = []
count = 0
for bed_file in bed_iterable:
    # clean name
    rep_w_blank = ".*/|.merged.sorted|.sorted|.bed$|.bed.gz$|.txt$|_gene_ct"
    bed_name = re.sub(rep_w_blank, "", bed_file)
    bed_name = re.sub("all_predictions", "cvdc_enhancers_dickel", bed_name)
    # use as column name
    print(bed_name)
    bed_df = pd.read_table(bed_file, names=bed_header)
    # print(bed_header)
    bed_df = bed_df[cols_to_keep]
    bed_df["anno"] = bed_name
    bed_df_list.append(bed_df)
    count += 1
    # if count > 3:
    #     break
    print(count)


bed_df = pd.concat(bed_df_list)

bed_df.columns
bed_df.head()
bed_df.shape

# saving
bed_df.to_csv(re.sub(".bed", "_anno_cts_grouped.txt", gene_bed_f),
              sep="\t")  # keep index after groupby index=False,
# _anno_cts.txt _anno_cts_grouped.txt
# mp result shape: (26624244, 1)
# non-mp results shape: (26624242, 1)

# sum over evertyhing
bed_df = bed_df.groupby(["Gene", "anno"]).sum()

# clean column names
# rep_w_blank = ".*/|.merged.sorted|.sorted|.bed$|.bed.gz$|.txt$"
# anno_list = [re.sub(rep_w_blank, "", i) for i in file_loc_list]
# anno_list = [re.sub("all_predictions", "cvdc_enhancers_dickel", i)
#              for i in anno_list]
