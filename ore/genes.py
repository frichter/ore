#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Gene processing.

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-01-22
:Copyright: 2018, Felix Richter
:License: CC BY-SA
"""

import os
import re
import subprocess

import pandas as pd

import pysam


class Error(Exception):
    """Base class for exceptions in this module."""

    pass


class RNASeqError(Error):
    """Exception raised for errors in this module.

    Attributes:
        message -- explanation of the error

    """

    def __init__(self, message):
        """Assign error explanation to object."""
        self.message = message


class Genes(object):
    """Object containing methods relevant to genes."""

    def __init__(self, gene_pheno_loc, var_bed_loc):
        """Initialize object containing methods relevant to genes.

        TODO:
            Check if file is tabix

        """
        if not os.path.exists(gene_pheno_loc):
            raise RNASeqError("This file does not exist: '{}'"
                              .format(gene_pheno_loc))
        self.gene_pheno_loc = gene_pheno_loc
        self.var_bed_loc = var_bed_loc
        tbx_handle = pysam.TabixFile(gene_pheno_loc)
        self.contigs = tbx_handle.contigs
        tbx_handle.close()
        self.check_gene_ref_genome()

    def assign_genes(self, current_chrom, upstream_only=False,
                     downstream_only=False, max_tss_dist=1e4,
                     gene_strand_data=None):
        """Assign variants to genes.

        Top level function for assigning genes
        variant input file prepared in vcf.py

        Assign variants to closest gene. only keep variants within X
            bp of gene TSS and update variant bed file as "nearTSS.bed"

        Docs for bedtools closest:
            http://bedtools.readthedocs.io/en/latest/content/tools/closest.html
            - D b option: reports distance and declare variant as
                upstream/downstream with respect to the strand
                of the gene TSS. Note that overlapping feature distance = 0
            -t option: what to do with genes with the same TSS? pick the first

        Args:
            upstream_only (:obj:`logical`): only variants UPstream of TSS
            downstream_only (:obj:`logical`): only variants DOWNstream of TSS
            max_tss_dist (:obj:`int`): maximum distance from TSS (used to
                increase efficiency by limiting size of data)

        TODO
            Is there a bedtools temp directory?

        """
        self.current_chrom = current_chrom  # "chr" +
        # create input files and declare file names
        gene_bed_loc = self.create_gene_per_chrom_file(gene_strand_data)
        var_bed_loc = self.var_bed_loc % self.current_chrom
        self.closest_gene_var_loc = re.sub("bed$", "closest_gene.bed",
                                           var_bed_loc)
        self.nearTSS_loc = re.sub("bed$", "nearTSS.bed", var_bed_loc)
        # Identify closest gene with bedtools
        bt_cmd = "time bedtools closest -a {} -b {} -D b -t first > {}".format(
            var_bed_loc, gene_bed_loc, self.closest_gene_var_loc)
        if not os.path.exists(self.closest_gene_var_loc):
            print(bt_cmd)
            subprocess.call(bt_cmd, shell=True)
        if not os.path.exists(self.nearTSS_loc):
            self.subset_variants_near_TSS(upstream_only, downstream_only,
                                          max_tss_dist)
        else:
            current_chrom = "Not_rerun_" + current_chrom
        return current_chrom

    def create_gene_per_chrom_file(self, gene_strand_data):
        """Create bed file containing gene location and strand.

        Load the gene strand information as `strand_dict` from
            `GENE_DF_LOCATION`, necessary because the FastQTL input file
            does not have strand information but this is important for
            determining upstream/downstream.

        Returns:
            `gene_bed_loc`: location of each gene's TSS along with strand

        TODO:
            Un-hardcode as much of this as possible (does ANNOVAR give strand?)
                e.g., GENE_DF_LOCATION
            remove genes with duplicate TSS (when preparing RNAseq)
            Why is `line_count` necessary?

        """
        gene_df = pd.read_table(gene_strand_data)
        strand_dict = dict(zip(gene_df.gene, gene_df.Strand))
        tbx_gene_handle = pysam.TabixFile(self.gene_pheno_loc)
        gene_bed_loc = re.sub(
            "tmp", "tmp_gene", self.var_bed_loc) % self.current_chrom
        line_count = 0
        with open(gene_bed_loc, 'w') as gene_bed_f:
            search_chrom = self.current_chrom
            if not self.ucsc_ref_genome:
                search_chrom = re.sub('chr', '', search_chrom)
            for line in tbx_gene_handle.fetch(search_chrom, 0, 3e8):
                line_list = line.strip().split("\t")[0:4]
                line_list.append(str(line_count))
                if line_list[3] in strand_dict:
                    line_list.append(strand_dict[line_list[3]])
                else:
                    line_list.append("NA")
                out_line = "\t".join(line_list)
                # if not self.ucsc_ref_genome:
                #     out_line = "chr" + out_line
                gene_bed_f.write(out_line + "\n")
                line_count += 1
        return gene_bed_loc

    def check_gene_ref_genome(self):
        """Determine reference genome used in gene phenotype file."""
        self.ucsc_ref_genome = False
        chrom_df = pd.read_table(self.gene_pheno_loc, nrows=5, usecols=[0])
        if chrom_df.iloc[0].astype(str).str.contains('chr').all():
            self.ucsc_ref_genome = True

    def subset_variants_near_TSS(self, upstream_only, downstream_only,
                                 max_tss_dist):
        """Take the subset of variants near the TSS.

        Only keep a subset of the variants, those within 10^5 base pairs of
            the TSS and possibly only those upstream or downstream.
        Can go under genes.py or annotations.py
        Note that filter is <= max_tss_dist and not < max_tss_dist so that a
            tss_distance cut-off of 0 corresponds to overlapping features

        Args:
            upstream_only: logical
            downstream_only: logical
            max_tss_dist: integer

        """
        cols_to_use = [0, 1, 2, 3, 4, 5, 8, 9, 11, 12]
        bed_col_names = ["Chrom", "Start", "End", "Ref", "Alt", "VCF_af",
                         "gene_TSS", "gene", "gene_strand", "tss_dist"]
        var_df = pd.read_table(self.closest_gene_var_loc,
                               header=None,
                               usecols=cols_to_use,
                               names=bed_col_names)
        # max TSS distance
        var_df = var_df.loc[abs(var_df.tss_dist) <= max_tss_dist]
        if upstream_only:
            var_df = var_df.loc[var_df.tss_dist < 0]
        elif downstream_only:
            var_df = var_df.loc[var_df.tss_dist > 0]
        # var_id should use 1-based start position
        var_df.Start1b = var_df.Start + 1
        var_df["var_id"] = (var_df.Chrom.str.cat(var_df.Start1b.
                            astype(str), sep='.').str.cat(var_df.Ref, sep='.').
                            str.cat(var_df.Alt, sep='.'))
        var_df.to_csv(self.nearTSS_loc, sep="\t", header=False, index=False,
                      float_format='%g')


# prepare gene file with TSS distances

#
