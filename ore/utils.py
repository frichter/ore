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
        raise ValueError


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
    # chrom_iter = itertools.chain([str(i) for i in range(1, 23)], ["X"])
    chrom_iter = itertools.chain([str(i) for i in contigs])
    # , "Y"
    pool = mp.Pool(processes=n_processes)
    # print("Total available cores: " + str(mp.cpu_count()))
    chroms_completed = pool.map(mp_function, chrom_iter)
    return chroms_completed


def applyParallel(dfGrouped, func):
    """Parallelize the pandas apply function.

    Source:
        https://stackoverflow.com/questions/26187759/parallelize-apply-after-pandas-groupby

    Args:
        dfGrouped (:obj:`DataFrame`): grouped dataframe, where
            `func` is to be applied to each group separately
        func (:obj:`int`):

    Returns:
        pd.concat(ret_list) (:obj:`DataFrame`): `dfGrouped` with function
            applied to every group

    """
    # from multiprocessing import Pool, cpu_count
    print("Using all {} cores".format(mp.cpu_count()))
    with mp.Pool(mp.cpu_count()) as p:
        ret_list = p.map(func, [group for name, group in dfGrouped])
    try:
        return pd.concat(ret_list)
    except ValueError:
        print("Was the dataframe being passed empty?", dfGrouped.head())
        raise ValueError


def anno_file_locations():
    """Import hardcoded annotation locations.

    Returns:
        file_loc_list (:obj:`list`): list of BED file annotations

    """
    data_dir = resource_filename('ore', 'data/hg19_genome_masks/')
    # mappability annotations
    rmsk = data_dir + "rmsk/rmsk.merged.sorted.bed"
    segdup = data_dir + "hg19_segdup.bed.gz"
    lcr = data_dir + "hg19_lcr_hs37d5.bed.gz"
    map300 = data_dir + "mappability300/mappability1_300.txt"
    hla_muc = data_dir + "genome/genes.MUC.HLA.bed"
    dac_blacklist = data_dir + "dac_blacklist.bed"
    duke_blacklist = data_dir + "encode_duke_blacklist.bed"
    pseudoauto_XY = data_dir + "pseudoautosomal_XY.bed"
    file_loc_list = [rmsk, segdup, lcr, map300, hla_muc, dac_blacklist,
                     duke_blacklist, pseudoauto_XY]
    file_loc_list = [segdup, lcr]
    return file_loc_list
