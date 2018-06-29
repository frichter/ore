#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Overlap variants with annotations.

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-01-22
:Copyright: 2018, Felix Richter
:License: CC BY-SA
"""


import os
import re
import subprocess
import glob

import pandas as pd

import pybedtools
from pybedtools import BedTool

from .utils import anno_file_locations


class Error(Exception):
    """Base class for exceptions in this module."""

    pass


class AnnotationError(Error):
    """Exception raised for errors in this module.

    Attributes:
        message (:obj:`str`): explanation of the error

    """

    def __init__(self, message):
        """Assign error explanation to object."""
        self.message = message


class Annotations(object):
    """Handles all annotation-related features.

    Attributes:
        use_annovar (:obj:`boolean`): T/F for whether or not ANNOVAR should
            be used to obtain allele frequencies and variant
            functional classes
        annovar_file_loc (:obj:`str`): ANNOVAR outputs with allele freqs
        var_bed_loc (:obj:`str`): variants in bed file
        long012_loc (:obj:`str`): variants in long format
        nearTSS_loc (:obj:`str`): variants within the maximum specified
            distance of any gene's TSS
        anno_out_loc (:obj:`str`): variants overlapped with annotations
        anno_out_list_loc (:obj:`str`): list of annotations that were
            overlapped with the variants
        final_var_loc (:obj:`str`): final set of variants (in 0/1/2 and long
            format) that will be used for association study, which are
            annotated with intra-cohort and GnomAD allele frequencies,
            variant functional class, and the closest gene
        annovar_dir (:obj:`str`): annovar command location
        humandb_dir (:obj:`str`): annovar database location
        genome_v (:obj:`str`): genome version for annovar (currently
            only hg19 or hg38)

    """

    def __init__(self, use_annovar, annovar_file_loc, var_bed_loc, long012_loc,
                 annovar_dir, humandb_dir, genome_v):
        """Create object containing methods/attributes for annotations.

        Args:
            use_annovar (:obj:`boolean`): T/F for whether or not ANNOVAR should
                be used to obtain allele frequencies and variant
                functional classes
            annovar_file_loc (:obj:`str`): output file with allele freqs
            var_bed_loc (:obj:`str`): variants in bed file
            long012_loc (:obj:`str`): variants in long format
            annovar_dir (:obj:`str`): annovar command location
            humandb_dir (:obj:`str`): annovar database location
            genome_v (:obj:`str`): genome version for annovar (currently
                only accepting hg19 or hg38, with automatic conversion
                of b37/hGRC38 to UCSC formats)

        """
        # intermediate file locations
        self.long012_loc, self.annovar_file_loc = long012_loc, annovar_file_loc
        self.nearTSS_loc = re.sub("bed$", "nearTSS.bed", var_bed_loc)
        self.anno_out_loc = re.sub("bed$", "anno.bed", var_bed_loc)
        self.anno_out_list_loc = re.sub(".bed$", "_list.txt",
                                        self.anno_out_loc)
        self.final_var_loc = re.sub("tmp_long_012", "var",
                                    self.long012_loc)
        # annovar info
        self.use_annovar = use_annovar
        self.annovar_dir, self.humandb_dir, self.genome_v = (
            annovar_dir, humandb_dir, genome_v)

    def run_annovar_cmd(self, current_chrom):
        """Execute ANNOVAR command.

        Input file prepared in vcf.py. Checks if annovar command has
            previously been run, checks if annovar should NOT be
            run, specifies the ANNOVAR command for 3 databases
            (refGene, ensGene, and gnomad_genome), and uses subprocess
            to execute the annovar command in the shell.

        Args:
            current_chrom (:obj:`str`): region of genome completed (b37 format)

        Returns:
            current_chrom (:obj:`str`): region of genome completed (b37 format)

        Raises:
            :obj:`AnnotationError`: if ANNOVAR command fails to execute

        """
        infile = self.annovar_file_loc % (current_chrom)  # "chr" +
        annovar_out_loc = re.sub(".txt", "", infile)
        if os.path.exists(annovar_out_loc + ".hg19_multianno.txt"):
            # print("ANNOVAR already done for %s" % infile)
            return "Not_rerun_" + current_chrom
        if not self.use_annovar:
            # print("ANNOVAR already done for %s" % infile)
            return "Not_run_" + current_chrom
        annovar_cmd = ("time perl {}table_annovar.pl " +
                       "{} {} -buildver {} " +
                       "-out {} --otherinfo -remove -protocol " +
                       "refGene,ensGene,gnomad_genome " +
                       # kaviar_20150923
                       "-operation g,g,f -nastring NA").format(
                       self.annovar_dir, infile,  self.humandb_dir,
                       self.genome_v, annovar_out_loc)
        print(annovar_cmd)
        subprocess.call(annovar_cmd, shell=True)
        return current_chrom

    def overlap_w_annotations(self, current_chrom):
        """Overlap variants with annotations.

        Args:
            current_chrom (:obj:`str`): region of genome for overlap

        Returns:
            current_chrom (:obj:`str`): region of genome completed (b37 format)

        Raises:
            :obj:`bed_error`: raised if bedtools errors out (most likely
                because of too little disk space for intermediate files)

        """
        # current_chrom = "chr" + current_chrom
        anno_out_loc = self.anno_out_loc % current_chrom
        anno_out_list_loc = self.anno_out_list_loc % current_chrom
        if os.path.exists(anno_out_loc):
            # print("Overlap w annotations already done for", current_chrom)
            return "Not_rerun_" + current_chrom
        file_loc_list = anno_file_locations()
        self.confirm_annotation_locations(file_loc_list, current_chrom)
        count = 0
        pybedtools.set_tempdir('pybedtools_temp_dir/')
        current_temp_dir = pybedtools.get_tempdir()
        var_bed_obj = BedTool(self.nearTSS_loc % current_chrom)
        overlap_list = []
        for file_loc in file_loc_list:
            bed_iterable = glob.iglob(file_loc)
            for bed_name in bed_iterable:
                overlap_list.append(bed_name)
                bed = BedTool(bed_name)
                try:
                    var_bed_obj = var_bed_obj.intersect(bed, c=True)
                except pybedtools.helpers.BEDToolsError as bed_error:
                    print("Check that the temp directory ({}) used by " +
                          "bedtools has space".format(current_temp_dir))
                    raise bed_error
                count += 1
                if count % 5 == 0:
                    print("Done with", count, bed_name)
        print("Completed overlaps for {}/{} annotations for {}".format(
            str(count), len(overlap_list), current_chrom))
        var_bed_obj.saveas(anno_out_loc)
        with open(anno_out_list_loc, 'w') as f:
            f.write("\n".join(overlap_list) + "\n")
        pybedtools.cleanup()
        return current_chrom

    def confirm_annotation_locations(self, file_loc_list, current_chrom):
        """Confirm that all files exist, count the number of files.

        Args:
            file_loc_list (:obj:`list`): list of files to overlap with variants
            current_chrom (:obj:`str`): region of genome for overlap

        Raises:
            :obj:`AnnotationError`: if the annotation files do not exist
            :obj:`bed_error`: raised if bedtools errors out while loading
                the annotation (most likely the annotation is not in correct
                bed format)

        """
        count = 0
        for file_loc in file_loc_list:
            file_iterable = glob.iglob(file_loc)
            for file_name in file_iterable:
                if os.path.isfile(file_name):
                    count += 1
                else:
                    raise AnnotationError("{} does not exist!".
                                          format(file_name))
                try:
                    BedTool(file_name)
                except pybedtools.helpers.BEDToolsError as bed_error:
                    print(
                        "Check that {} is in correct BED format".format(
                            file_name))
                    raise bed_error
        print(
            "Overlapping {} annotation files with variants on {}".format(
                count, current_chrom))

    def get_final_set_of_variants(self, current_chrom):
        """Process ANNOVAR, annotation overlap, and 012 results.

        Args:
            current_chrom (:obj:`str`): region of genome being considered

        Returns:
            current_chrom (:obj:`str`): region of genome completed (b37 format)

        """
        # current_chrom = "chr" + current_chrom
        if os.path.exists(self.final_var_loc % current_chrom):
            # print("Final variants already prepared for %s" % current_chrom)
            return "Not_rerun_" + current_chrom
        print("Cleaning Annotations for", current_chrom)
        anno_df = self.import_vars_close_to_gene(current_chrom)
        # print("Closest gene DF size:", anno_df.shape)
        anno_df = anno_df.set_index('var_id')
        if self.use_annovar:
            print("Cleaning ANNOVAR results for", current_chrom)
            annovar_df = self.clean_annovar_results(current_chrom)
            # print("ANNOVAR DF size:", annovar_df.shape)
            print("Joining ANNOVAR with annotations for", current_chrom)
            anno_df = annovar_df.set_index('var_id').join(anno_df, how='inner')
            # print("ANNOVAR joined w closest gene DF size:", anno_df.shape)
        clean_df = self.remove_vars_in_unwanted_cols(anno_df)
        # print("Repeats and segdups removed DF size:", clean_df.shape)
        print("Loading long 012 matrix for", current_chrom)
        long012_df = self.load_long_012_df(current_chrom)
        # print("012 long DF size:", long012_df.shape)
        print("Joining long 012 with annotations", current_chrom)
        final_df = clean_df.join(long012_df.set_index('var_id'), how='inner')
        # print("012 long joined with annotated DF size:", final_df.shape)
        print("Getting intra-cohort variant counts/frequency")
        final_df.reset_index(inplace=True)
        final_df['var_id_count'] = final_df.groupby(
            'var_id')['GT'].transform('sum')
        # multiply by 2 for autosomes, 1 for sex chromosomes in MALES
        # This should be reflected in the maximum genotype that can
        # be observed for a given chromosome.. however there seem
        # to be plenty of Y chromosome variants with GT=1 and GT=2...
        id_ct = len(set(final_df.blinded_id)) * 2
        final_df['var_id_freq'] = final_df.var_id_count/id_ct
        if not self.use_annovar:
            print("Setting popmax AF to 0 (because not using ANNOVAR)")
            final_df['popmax_af'] = 0
        print("Writing to", self.final_var_loc % current_chrom)
        final_df.to_csv(self.final_var_loc % current_chrom, sep="\t",
                        index=False, float_format='%g')
        print("Done writing to", self.final_var_loc % current_chrom)
        return current_chrom

    def clean_annovar_results(self, current_chrom):
        """Process ANNOVAR output.

        calculate pop_max_1k (using populations with n>1k on GnomAD) and
            pop_max (using all populations) and decide which other
            ANNOVAR info should be joined with 012 file (gene,
            var_function, tss_distance) (ie ignore ChooseBestGene
            function in rare_variant_prep.R)
        Create bed file of variants with pop_max_1k <0.01 (for eQTLs)
        Join ANNOVAR output with 012 files (here or in variants.py)

        TODO:
            confirm annovar_df.shape rows is >= anno_df.shape rows and that
            rows in anno_df.shape is == joined_anno_df.shape
            if the latter is not true, run this command:
            anno_df[~anno_df.var_id.isin(annovar_df.var_id)]

        Args:
            current_chrom (:obj:`str`): region of genome being considered

        Returns:
            current_chrom (:obj:`str`): region of genome completed (b37 format)

        TODO:
            low_memory is deprecated so instead specify any ambiguous dtypes
            https://stackoverflow.com/questions/24251219/pandas-read-csv-low-memory-and-dtype-options

        """
        annovar_out_loc = re.sub("txt$", "hg19_multianno.txt",
                                 self.annovar_file_loc % current_chrom)
        try:
            annovar_df = pd.read_table(annovar_out_loc, low_memory=False)
        except pd.errors.ParserError as pandas_import_error:
                print("")
                raise pandas_import_error
        # rename and keep only some of the columns
        col_rename_dict = {"Otherinfo": "var_id",
                           "Func.refGene": "func_refgene",
                           "Gene.refGene": "gene_refgene",
                           "GeneDetail.refGene": "dist_refgene",
                           "ExonicFunc.refGene": "exon_func_refgene",
                           "Func.ensGene": "func_ensgene",
                           "Gene.ensGene": "gene_ensgene",
                           "GeneDetail.ensGene": "dist_ensgene",
                           "ExonicFunc.ensGene": "exon_func_ensgene"}
        annovar_df = annovar_df.rename(columns=col_rename_dict)
        cols_to_keep = ["var_id", "func_refgene", "gene_refgene",
                        "dist_refgene", "exon_func_refgene",
                        "func_ensgene", "gene_ensgene", "dist_ensgene",
                        "exon_func_ensgene",
                        "gnomAD_genome_ALL",
                        "gnomAD_genome_AFR", "gnomAD_genome_FIN",
                        "gnomAD_genome_NFE"]
        annovar_df = annovar_df.reindex(columns=cols_to_keep)
        # popmax using only populations with n>1000
        popmax_col_list = ["gnomAD_genome_ALL", "gnomAD_genome_AFR",
                           "gnomAD_genome_FIN", "gnomAD_genome_NFE"]
        annovar_df.loc[:, popmax_col_list] = (
            annovar_df.loc[:, popmax_col_list].
            apply(pd.to_numeric, errors='coerce').fillna(0))
        if pd.isnull(annovar_df.gnomAD_genome_ALL).any():
            raise TypeError("Unable to convert NaNs to 0 in AF columns")
        annovar_df["popmax_af"] = annovar_df[popmax_col_list].max(axis=1)
        annovar_df.drop(popmax_col_list, axis=1, inplace=True)
        return annovar_df

    def import_vars_close_to_gene(self, current_chrom):
        """Import variants that are close to the TSS.

        Import and clean column names, then import the actual bed file.

        Args:
            current_chrom (:obj:`str`): region of genome being considered

        Returns:
            current_chrom (:obj:`str`): region of genome completed (b37 format)

        TODO:
            Generalize filter to other annotations (i.e., not just enhancers)

        """
        anno_out_list_loc = self.anno_out_list_loc % current_chrom
        anno_out_loc = self.anno_out_loc % current_chrom
        anno_list = pd.read_table(anno_out_list_loc, squeeze=True, header=None)
        # clean column names
        rep_w_blank = ".*/|.merged.sorted|.sorted|.bed$|.bed.gz$|.txt$"
        anno_list = [re.sub(rep_w_blank, "", i) for i in anno_list]
        anno_col_names = ["Chrom", "Start0", "End0", "Ref", "Alt", "VCF_af",
                          "gene_TSS", "gene", "gene_strand", "tss_dist",
                          "var_id"]
        anno_col_names.extend(anno_list)
        anno_df = pd.read_table(anno_out_loc, header=None,
                                names=anno_col_names)
        # filter for only enhancers (for now)
        # anno_df = anno_df[anno_df.cvdc_enhancers_dickel > 0]
        return anno_df

    def remove_vars_in_unwanted_cols(self, joined_anno_df):
        """Remove variants in repeats and segdups.

        Args:
            joined_anno_df (:obj:`DataFrame`): pandas dataframe with
                column names that include those listed in
                unwanted_cols

        Raises:
            :obj:`AnnotationError`: if variants have not been overlapped with
                regions that you want to exclude

        TODO:
            Un-hardcode this/customize which columns to remove

        """
        unwanted_cols = ['hg19_segdup', 'hg19_lcr_hs37d5']
        try:
            unwanted_vars_df = joined_anno_df.loc[:, unwanted_cols] == 0
        except KeyError:
            raise AnnotationError(
                "In order to remove variants in regions, please first " +
                "overlap the variants {}".format(", ".join(unwanted_cols)))
        else:
            clean_df = joined_anno_df[unwanted_vars_df.all(axis=1)]
        return clean_df

    def load_long_012_df(self, current_chrom):
        """Load matrix containing genotype status per variant per person.

        Args:
            current_chrom (:obj:`str`): region of genome being considered

        Returns:
            long012_df (:obj:`DataFrame`): pandas dataframe with 3 columns,
                blinded_id, GT, and var_id (the latter being a joined string
                of chrom.pos.ref.alt)

        """
        long_col_names = ['loc_id', 'blinded_id', 'GT', 'Ref', 'Alt']
        long012_df = pd.read_table(self.long012_loc % current_chrom,
                                   header=None, names=long_col_names)
        long012_df["var_id"] = (long012_df.loc_id.str.
                                cat(long012_df.Ref.astype(str), sep='.').
                                str.cat(long012_df.Alt, sep='.'))
        long012_df.drop(['Ref', 'Alt', 'loc_id'], axis=1, inplace=True)
        # long012_df = long012_df[long012_df.GT == 1]
        return long012_df
