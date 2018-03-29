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
    # logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    return(logger)


def prepare_directory(new_dir, clean_run=False):
    """Prepare per chromosome directory FOR NEW RUN.

    Args:
        new_dir (:obj:`str`): Location of new directory to create

    Create a per_chrom directory if it does not exist. If it does exist,
        ask the user if the directory should be cleaned for
        a new run.

    """
    print("Working in this directory:", os.getcwd())
    if not os.path.exists(new_dir):
        print("Creating", new_dir)
        os.makedirs(new_dir)
    if clean_run:
        [os.remove(f) for f in glob.iglob(new_dir + "/*")]


def applyParallel(dfGrouped, func):
    """Parallelize the pandas apply function.

    SOURCE
        https://stackoverflow.com/questions/26187759/parallelize-apply-after-pandas-groupby

    TODO
        handle ValueError: No objects to concatenate
    """
    from multiprocessing import Pool, cpu_count
    print("Using all {} cores".format(cpu_count()))
    with Pool(cpu_count()) as p:
        ret_list = p.map(func, [group for name, group in dfGrouped])
    return pd.concat(ret_list)


def anno_file_locations():
    """Import hardcoded annotation locations."""
    data_dir = __file__ + "/data/"
    # mappability annotations
    rmsk = data_dir + "rmsk/rmsk.merged.sorted.bed"
    segdup = data_dir + "segdup/segdup.merged.sorted.bed"
    lcr = data_dir + "LCR-hs37d5_chr.bed"
    map300 = data_dir + "mappability300/mappability1_300.txt"
    hla_muc = data_dir + "genome/genes.MUC.HLA.bed"
    dac_blacklist = data_dir + "dac_blacklist.bed"
    duke_blacklist = data_dir + "encode_duke_blacklist.bed"
    pseudoauto_XY = data_dir + "pseudoautosomal_XY.bed"
    file_loc_list = [rmsk, segdup, lcr, map300, hla_muc, dac_blacklist,
                     duke_blacklist, pseudoauto_XY]
    return file_loc_list
