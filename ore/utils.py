#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Various utilities.

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-03-26
:Copyright: 2018, Felix Richter
:License: CC BY-SA
"""

import logging
import sys
import os
import glob
import multiprocessing as mp
import itertools

from pkg_resources import resource_filename
import pandas as pd


def initialize_logger(log_file, logAppName):
    """Set up logging.

    Returns:
        logger (:obj:`logging object`): log with file and stream handlers

    """
    logger = logging.getLogger(logAppName)
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG)
    # create formatter and add it to the handlers
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(formatter)
    formatter = logging.Formatter(
        "%(asctime)s - %(message)s")
    console_handler.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    return(logger)


def checkCPUcount(n_processes):
    """Confirm CPU count input by user is less than max cores and int.

    Args:
        n_processes (:obj:`int`): number of cores that the user wants to use

    """
    if n_processes > mp.cpu_count():
        print("Current machine only has {} cores, but input was {}".format(
              mp.cpu_count(), n_processes))
        # raise ValueError
        n_processes = mp.cpu_count() - 1
    return n_processes


def check_variant_inputs(args):
    """Confirm variant-related inputs.

    Args:
        args (:obj:`argparse object`): contains user inputs

    """
    # TSS distance > 0
    if any([i <= 0 for i in args.tss_dist]):
        print("TSS distance must be greater than 0")
        raise ValueError
    if args.upstream and args.downstream:
        print("Cannot use BOTH --upstream and --downstream")
        raise ValueError


def check_ANNOVAR_inputs(args):
    """Confirm that ANNOVAR related inputs are available."""
    if args.annovar:
        if not os.path.exists(args.annovar_dir + "table_annovar.pl"):
            print("table_annovar.pl does not exist in", args.annovar_dir)
            raise FileNotFoundError
        if not os.path.exists(args.annovar_dir + "annotate_variation.pl"):
            print("annotate_variation.pl does not exist in", args.annovar_dir)
            raise FileNotFoundError
        # list of inputs needed for ANNOVAR in humandb directory
        # provide code to install if not available
        # what inputs are needed in humandb?


def prepare_directory(new_dir, clean_run=False):
    """Prepare per chromosome directory FOR NEW RUN.

    Args:
        new_dir (:obj:`str`): Location of new directory to create

    Create a per_chrom directory if it does not exist. If it does exist,
        ask the user if the directory should be cleaned for
        a new run.

    """
    # print("Working in this directory:", os.getcwd())
    if not os.path.exists(new_dir):
        print("Creating", new_dir)
        os.makedirs(new_dir)
    if clean_run:
        [os.remove(f) for f in glob.iglob(new_dir + "/*")]


def multiprocess_by_chrom_cmd(n_processes, contigs, mp_function):
    """Loop over chromosomes (by sending to multiple processes).

    Args:
        n_processes (:obj:`int`): number of workers/cores to run at a time
        contigs (:obj:`list`): contigs present in the VCF
        mp_function (:obj:`function`): function being processed across
            multiple genomic regions

    Returns:
        chroms_completed (:obj:`list`): genomic regions that were
            processed by `mp_function`

    """
    chrom_iter = itertools.chain([str(i) for i in contigs])
    pool = mp.Pool(processes=n_processes)
    print("Using {} out of {} cores".format(n_processes, mp.cpu_count()))
    chroms_completed = pool.map(mp_function, chrom_iter)
    # chroms_completed = []
    # for chrom_i in chrom_iter:
    #     chrom_done = mp_function(chrom_i)
    #     chroms_completed.append(chrom_done)
    return chroms_completed


def applyParallel(dfGrouped, func, n_processes):
    """Parallelize the pandas apply function.

    Source:
        https://stackoverflow.com/questions/26187759/parallelize-apply-after-pandas-groupby

    Args:
        dfGrouped (:obj:`DataFrame`): grouped dataframe, where
            `func` is to be applied to each group separately
        func (:obj:`int`):
        n_processes (:obj:`int`): number of workers/cores to run at a time

    Returns:
        pd.concat(ret_list) (:obj:`DataFrame`): `dfGrouped` with function
            applied to every group

    """
    # from multiprocessing import Pool, cpu_count
    print("Using {} out of {} cores".format(n_processes, mp.cpu_count()))
    # print("Total available cores: " + str(mp.cpu_count()))
    with mp.Pool(processes=n_processes) as p:
        ret_list = p.map(func, [group for name, group in dfGrouped])
    try:
        return pd.concat(ret_list)
    except ValueError:
        print("Was the dataframe being passed empty?", dfGrouped.head())
        raise ValueError


def anno_file_locations(annotations):
    """Import hardcoded annotation locations.

    Returns:
        file_loc_list (:obj:`list`): list of BED file annotations

    """
    data_dir = resource_filename('ore', 'data/hg19_genome_masks/')
    # mappability annotations
    # rmsk = data_dir + "rmsk/rmsk.merged.sorted.bed"
    segdup = data_dir + "hg19_segdup.bed.gz"
    lcr = data_dir + "hg19_lcr_hs37d5.bed.gz"
    # map300 = data_dir + "mappability300/mappability1_300.txt"
    # hla_muc = data_dir + "genome/genes.MUC.HLA.bed"
    # dac_blacklist = data_dir + "dac_blacklist.bed"
    # duke_blacklist = data_dir + "encode_duke_blacklist.bed"
    # pseudoauto_XY = data_dir + "pseudoautosomal_XY.bed"
    # anno_hg19 = ("/hpc/users/richtf01/chdiTrios/Felix/wgs/bed_annotations/" +
    #              "hg19_all/*.bed")
    # file_loc_list = [rmsk, segdup, lcr, map300, dac_blacklist,  # hla_muc,
    #                  duke_blacklist, pseudoauto_XY]
    file_loc_list = [segdup, lcr]
    # file_loc_list = [anno_hg19]
    # file_loc_list = [
    #     "hg19_mapping_2018_07/*bed",
    #     "deepheart_anno/human_data/dickel_2016_hg19/*bed",
    #     "deepheart_anno/human_data/dickel_2016_lifted/*lifted_hg19.bed",
    #     "deepheart_anno/human_data/encode/*.bed",
    #     "deepheart_anno/human_data/roadmap_epigenomics/*.narrowPeak",
    #     "deepheart_anno/mouse_lifted_hg19/*.bed",
    #     "deepheart_anno/human_data/uw_gsms/*.bed"]
    # file_loc_list = ["/hpc/users/richtf01/chdiTrios/Felix/wgs/" +
    #                  "bed_annotations/" + i for i in file_loc_list]
    if annotations:
        file_loc_list.extend(annotations)
    return file_loc_list


def filter_refgene_ensgene(var_df_per_chrom, variant_class,
                           refgene, ensgene):
    """Filter for a refgene function, ensembl function or both."""
    variant_class = "^" + variant_class
    if not refgene and not ensgene:
        print("Using RefGene for filtering")
        refgene = True
    if refgene:
        vars_refgene = var_df_per_chrom.func_refgene.str.contains(
            variant_class, regex=True)
        var_df_per_chrom = var_df_per_chrom[vars_refgene]
    if ensgene:
        vars_ensgene = var_df_per_chrom.func_ensgene.str.contains(
            variant_class, regex=True)
        var_df_per_chrom = var_df_per_chrom[vars_ensgene]
    return var_df_per_chrom


def filter_refgene_ensgene_exon(var_df_per_chrom, exon_class,
                                refgene, ensgene):
    """Filter for a refgene function, ensembl function or both.

    Args:
        var_df_per_chrom (:obj:`DataFrame`): all variants in a chromosome
        variant_class (:obj:`str`): annovar variant class to filter
            on (default None)
        exon_class (:obj:`str`): annovar EXON class to filter
            on (default None)
        refgene (:obj:`boolean`): if used RefSeq to define variant classes
        ensgene (:obj:`boolean`): using ENSEMBL to define variant classes

    Returns:
        var_df_per_chrom (:obj:`DataFrame`): only variants in the
            desired `exon_class`

    Description:
        First prepends a ^ so that only the highest impact `exon_class`
        is considered as the de-facto class for filtering.

    """
    exon_class = "^" + exon_class
    if not refgene and not ensgene:
        print("Using RefGene for filtering")
        refgene = True
    if refgene:
        vars_refgene = var_df_per_chrom.exon_func_refgene.str.contains(
            exon_class, regex=True)
        var_df_per_chrom = var_df_per_chrom[vars_refgene]
    if ensgene:
        vars_ensgene = var_df_per_chrom.exon_func_ensgene.str.contains(
            exon_class, regex=True)
        var_df_per_chrom = var_df_per_chrom[vars_ensgene]
    return var_df_per_chrom


def filter_variant_class(df, variant_class, exon_class, refgene, ensgene):
    """Filter a dataframe for the desired variant classes."""
    if variant_class:
        df = filter_refgene_ensgene(
            df, variant_class, refgene, ensgene)
        if variant_class.startswith("exonic") and exon_class:
            df = filter_refgene_ensgene_exon(
                df, exon_class, refgene, ensgene)
    return df
