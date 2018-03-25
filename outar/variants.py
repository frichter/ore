#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Process variants for enrichment.

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-01-20
:Copyright: 2018, Felix Richter
:License: CC BY-SA
"""

import os
import re

from functools import partial
import multiprocessing as mp
import itertools

from vcf import VCF
from annotations import Annotations
from genes import Genes
from outliers import Outliers


def prepare_per_chrom_directory(chrom_dir, vcf_obj):
    """Prepare per chromosome directory FOR NEW RUN.

    Create a per_chrom directory if it does not exist. If it does exist,
        remove previously generated per_chrom files (bc code appends)

    """
    print("Working in this directory:", os.getcwd())
    if not os.path.exists(chrom_dir):
        print("Creating", chrom_dir)
        os.makedirs(chrom_dir)
    os.chdir(chrom_dir)
    if not os.path.exists("pybedtools_temp_dir/"):
        print("Creating pybedtools_temp_dir/")
        os.makedirs("pybedtools_temp_dir/")


def multiprocess_by_chrom_cmd(n_processes, mp_function):
    """Loop over chromosomes (by sending to multiple processes).

    TODO:
        Non-human chromosomes (brainstorm how to obtain from file. tabix?)

    """
    chrom_iter = itertools.chain([str(i) for i in range(1, 23)], ["X"])
    # , "Y"
    # chrom_iter = itertools.chain([str(i) for i in range(1, 2)])
    pool = mp.Pool(processes=n_processes)
    # print("Total available cores: " + str(mp.cpu_count()))
    chroms_completed = pool.map(mp_function, chrom_iter)
    return chroms_completed


class Variants(object):
    """Instantiate object for processing Variants."""

    def __init__(self, vcf_loc, gene_pheno_loc, prefix, outlier_postfix=None,
                 annovar_dir=__file__ + '/annovar/',
                 humandb_dir=__file__ + '/annovar/humandb/',
                 n_processes=1):
        """Instantiate Variants class.

        Args:
            vcf_loc (:obj:`str`): VCF location
            n_processes (:obj:`int`): number of processes to use

        Attributes:
            vcf_obj: (:obj:`VCF`): VCF instance to maintain file locations

        """
        self.n_processes = n_processes
        self.vcf_loc = vcf_loc
        self.vcf_obj = VCF(self.vcf_loc, prefix=prefix)
        # make a common temporary directory
        self.per_chrom_dir = self.vcf_obj.prefix + "_per_chrom/"
        current_chrom_file_loc = (self.per_chrom_dir +
                                  "tmp_long_012_%s.txt")
        self.vcf_obj.declare_output_file_names(current_chrom_file_loc)
        print("VCF of interest loaded, with IDs parsed successfully")
        self.anno_obj = Annotations(self.vcf_obj.annovar_file_loc,
                                    self.vcf_obj.bed_file_loc,
                                    annovar_dir, humandb_dir,
                                    "hg19",
                                    current_chrom_file_loc)
        print("Annotation functions loaded...")
        self.gene_obj = Genes(gene_pheno_loc, self.vcf_obj.bed_file_loc)
        print("Gene object loaded...")
        self.outlier_obj = Outliers(gene_pheno_loc)
        print("Outliers initialized...")
        self.expr_outs_loc = (self.vcf_obj.prefix + "_outliers.txt")
        if outlier_postfix:
            self.expr_outs_loc = (self.vcf_obj.prefix + "_" + outlier_postfix)
        # self.expr_outs_loc = "wgs_singletons_n83_outs_norm_ext_rm1500.txt"

    def extract_variants_from_vcf(self):
        """Obtain allele-ID pairs from VCF in long format.

        Attributes:

        Returns:
            chroms_completed (:obj:`list`): chromosomes converted to
            long format (or failed)

        """
        # alternatively: check to see if files are complete (and don't delete)
        prepare_per_chrom_directory(self.per_chrom_dir, self.vcf_obj)
        # assign the constant variables to `prepare_vcf_per_chrom` in
        # preparation for multiprocessing
        partial_prepare_vcf_per_chrom = partial(
            self.vcf_obj.prepare_vcf_per_chrom, vcf_loc=self.vcf_loc,
            current_chrom_file_loc=self.vcf_obj.current_chrom_file_loc)
        chroms_completed = multiprocess_by_chrom_cmd(
            self.n_processes, partial_prepare_vcf_per_chrom)
        return chroms_completed

    def get_allele_frequencies(self):
        """Run the corresponding ANNOVAR command.

        Decide if you want run_annovar_function to be a static method

        """
        chroms_completed = multiprocess_by_chrom_cmd(
            self.n_processes, self.anno_obj.run_annovar)
        return chroms_completed

    def label_with_closest_gene(self, upstream_only, downstream_only,
                                max_tss_dist):
        """Find the closest gene."""
        partial_assign_genes = partial(
            self.gene_obj.assign_genes, upstream_only=upstream_only,
            downstream_only=downstream_only, max_tss_dist=max_tss_dist)
        chroms_completed = multiprocess_by_chrom_cmd(
            self.n_processes, partial_assign_genes)
        return chroms_completed

    def overlap_w_annotations(self):
        """Overlap variants with noncoding annotations."""
        chroms_completed = multiprocess_by_chrom_cmd(
            self.n_processes, self.anno_obj.overlap_vars_w_annotations)
        return chroms_completed

    def finalize_variants(self):
        """Get the final set of variants."""
        chroms_completed = multiprocess_by_chrom_cmd(
            self.n_processes, self.anno_obj.get_final_set_of_variants)
        return chroms_completed

    def prepare_outliers(self, most_extreme, distribution, outlier_max=None):
        """Obtain gene expression outliers."""
        if os.path.exists(self.expr_outs_loc):
            return "Already made outlier file " + self.expr_outs_loc
        # confirm IDs in WGS match RNAseq
        self.vcf_obj.load_vcf()
        ids_to_keep = self.vcf_obj.id_list
        expr_outlier_df = self.outlier_obj.get_outliers(
            ids_to_keep, distribution=distribution, most_extreme=most_extreme,
            expr_cut_off=2)
        outs_per_id_file = re.sub('.txt', '_outliers_per_id_ALL',
                                  self.expr_outs_loc)
        self.outlier_obj.plot_out_per_id(expr_outlier_df, outs_per_id_file)
        # determine which IDs have too many outliers (and remove these)
        if outlier_max:
            ids_to_keep = self.outlier_obj.get_ids_w_low_out_ct(
                expr_outlier_df, outlier_max=outlier_max)
            # call outliers with these IDs
            expr_outlier_df = self.outlier_obj.get_outliers(
                ids_to_keep, distribution=distribution,
                most_extreme=most_extreme, expr_cut_off=2)
            outs_per_id_file = re.sub('.txt', '_outliers_per_id_final',
                                      self.expr_outs_loc)
            self.outlier_obj.plot_out_per_id(expr_outlier_df, outs_per_id_file)
        # write expr_outlier_df to file
        print("Saving outlier status dataframe to", self.expr_outs_loc)
        expr_outlier_df.to_csv(self.expr_outs_loc, sep="\t", index=False)
