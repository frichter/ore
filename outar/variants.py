#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Process variants for enrichment.

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-01-20
:Copyright: 2018, Felix Richter
:License: CC BY-SA
"""

import os

from functools import partial

from .vcf import VCF
from .annotations import Annotations
from .genes import Genes
from .utils import prepare_directory, multiprocess_by_chrom_cmd


class Variants(object):
    """Instantiate object for processing Variants."""

    def __init__(self, vcf_loc, gene_pheno_loc,
                 output_prefix,
                 outlier_postfix,
                 use_annovar,
                 annovar_dir,
                 humandb_dir,
                 n_processes,
                 clean_run,
                 logger):
        """Instantiate Variants class.

        Args:
            vcf_loc (:obj:`str`): VCF location
            gene_pheno_loc (:obj:`str`): Expression file location
            output_prefix (:obj:`str`): file prefix for outputs
            outlier_postfix (:obj:`str`): file ending for outlier files
            annovar_dir (:obj:`str`): ANNOVAR table_annovar.pl script directory
            humandb_dir (:obj:`str`): ANNOVAR data directory
            n_processes (:obj:`int`): number of processes to use
            clean_run (:obj:`boolean`): if true, delete info from previous runs
            logger (:obj:`logging object`): Current logger

        Attributes:
            vcf_loc (:obj:`str`): VCF location
            vcf_obj: (:obj:`VCF`): VCF instance to maintain file locations
            anno_obj: (:obj:`Annotations`): class instance of all
                annotation relevant methods and states
            gene_obj: (:obj:`Genes`): class instance of methods and
                states relevant to gene positions
            n_processes (:obj:`int`): number of processes to use

        """
        self.n_processes = n_processes
        self.clean_run = clean_run
        self.vcf_loc = vcf_loc
        self.vcf_obj = VCF(self.vcf_loc, prefix=output_prefix)
        # make a common temporary directory
        self.per_chrom_dir = output_prefix + "_per_chrom/"
        current_chrom_file_loc = (self.per_chrom_dir +
                                  "tmp_long_012_%s.txt")
        self.vcf_obj.load_vcf()
        self.vcf_obj.declare_output_file_names(current_chrom_file_loc)
        logger.info("VCF of interest loaded, with IDs parsed successfully")
        self.anno_obj = Annotations(use_annovar,
                                    self.vcf_obj.annovar_file_loc,
                                    self.vcf_obj.bed_file_loc,
                                    current_chrom_file_loc,
                                    annovar_dir,
                                    humandb_dir,
                                    "hg19")
        logger.info("Annotation functions loaded...")
        self.gene_obj = Genes(gene_pheno_loc, self.vcf_obj.bed_file_loc)
        logger.info("Gene object loaded...")

    def extract_variants_from_vcf(self):
        """Obtain allele-ID pairs from VCF in long format.

        Attributes:

        Returns:
            chroms_completed (:obj:`list`): chromosomes converted to
                long format (or failed)

        """
        # alternatively: check to see if files are complete (and don't delete)
        prepare_directory(self.per_chrom_dir, self.clean_run)
        os.chdir(self.per_chrom_dir + "/..")
        prepare_directory("pybedtools_temp_dir/")
        # assign the constant variables to `prepare_vcf_per_chrom` in
        # preparation for multiprocessing
        partial_prepare_vcf_per_chrom = partial(
            self.vcf_obj.prepare_vcf_per_chrom,
            vcf_loc=self.vcf_loc,
            current_chrom_file_loc=self.vcf_obj.current_chrom_file_loc)
        chroms_completed = multiprocess_by_chrom_cmd(
            self.n_processes,
            self.vcf_obj.contigs,
            partial_prepare_vcf_per_chrom)
        return chroms_completed

    def run_annovar_wrapper(self):
        """Run the corresponding ANNOVAR command.

        Returns:
            chroms_completed (:obj:`list`): chromosomes converted to
                long format (or failed)

        """
        chroms_completed = multiprocess_by_chrom_cmd(
            self.n_processes,
            self.vcf_obj.contigs,
            self.anno_obj.run_annovar_cmd)
        return chroms_completed

    def label_with_closest_gene(self, upstream_only, downstream_only,
                                max_tss_dist, gene_strand_data):
        """Find the closest gene.

        Args:
            upstream_only (:obj:`bool`): If true, only considering variants
                upstream of the transcription start site (TSS)
            downstream_only (:obj:`bool`): If true, only considering variants
                DOWNstream of the TSS
            max_tss_dist (:obj:`int`): only consider variants within this
                distance of the TSS

        Returns:
            chroms_completed (:obj:`list`): chromosomes converted to
                long format (or failed)

        """
        partial_assign_genes = partial(
            self.gene_obj.assign_genes, upstream_only=upstream_only,
            downstream_only=downstream_only, max_tss_dist=max_tss_dist,
            gene_strand_data=gene_strand_data)
        chroms_completed = multiprocess_by_chrom_cmd(
            self.n_processes,
            self.vcf_obj.contigs,
            partial_assign_genes)
        return chroms_completed

    def overlap_w_annotations_wrapper(self):
        """Overlap variants with noncoding annotations.

        Returns:
            chroms_completed (:obj:`list`): chromosomes converted to
                long format (or failed)

        """
        chroms_completed = multiprocess_by_chrom_cmd(
            self.n_processes,
            self.vcf_obj.contigs,
            self.anno_obj.overlap_w_annotations)
        return chroms_completed

    def finalize_variants(self):
        """Get the final set of variants.

        Returns:
            chroms_completed (:obj:`list`): chromosomes converted to
                long format (or failed)

        """
        chroms_completed = multiprocess_by_chrom_cmd(
            self.n_processes,
            self.vcf_obj.contigs,
            self.anno_obj.get_final_set_of_variants)
        return chroms_completed
