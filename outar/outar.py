#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Prepare variants for enrichment.

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-01-20
:Copyright: 2018, Felix Richter
:License: CC BY-SA
"""


import argparse
import logging
import sys

from variants import Variants
from enrichment import Enrich

"""Profiling libraries:
import cProfile
import pstats
from memory_profiler import profile
"""


def initialize_logger(log_file):
    """Set up logging.

    Returns:
        logger (:obj:`logging object`): log with file and stream handlers

    """
    logger = logging.getLogger("OutAR_status")
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


# @profile
def associate_outliers(vcf_loc, gene_pheno_loc, output_prefix, outlier_postfix,
                       extrema, distribution, variant_class, enrich_loc):
    """Prepare and associate variants and outliers.

    Args:
        vcf_loc (:obj:`str`): VCF location
        gene_pheno_loc (:obj:`str`): gene expression (phenotype) location
        output_prefix (:obj:`str`): file prefix for outputs
        outlier_postfix (:obj:`str`): file prefix for outlier files
        extrema (:obj:`boolean`): T/F for using most extreme outlier
        distribution (:obj:`str`): type of outlier distribution considered
        variant_class (:obj:`str`): variant functional class to test
        enrich_loc (:obj:`str`): file location of enrichment results

    """
    logger = initialize_logger(log_file=output_prefix + "_outar_2.log")
    variants_obj = Variants(vcf_loc, gene_pheno_loc, prefix=output_prefix,
                            outlier_postfix=outlier_postfix, n_processes=12)
    chroms_completed = variants_obj.extract_variants_from_vcf()
    logger.info("LONG format variant abstraction done for...\n" +
                ", ".join(chroms_completed) + "\n")
    # annovar
    chroms_completed = variants_obj.get_allele_frequencies()
    logger.info("ANNOVAR done for...\n" +
                ", ".join(chroms_completed) + "\n")
    # find closest gene
    chroms_completed = variants_obj.label_with_closest_gene(
        upstream_only=False, downstream_only=False, max_tss_dist=1e4)
    logger.info("Closest gene found for the following chromosomes..\n" +
                ", ".join(chroms_completed) + "\n")
    # other annotations
    chroms_completed = variants_obj.overlap_w_annotations()
    logger.info("Overlaps with other bed files done for...\n" +
                ", ".join(chroms_completed) + "\n")
    # obtain outlier dataframe (write to file)
    variants_obj.prepare_outliers(most_extreme=extrema,
                                  distribution=distribution,
                                  outlier_max=1000)
    logger.info("Outliers prepared")
    # finalize variants (join variants with gene expression locations)
    chroms_completed = variants_obj.finalize_variants()
    logger.info("Final set of variants done for...\n" +
                ", ".join(chroms_completed) + "\n")
    # output final set of outliers and calculate enrichment
    rv_outlier_loc = output_prefix + "_rv_w_outliers.txt"
    enrich_obj = Enrich(variants_obj.anno_obj.final_var_loc,
                        variants_obj.expr_outs_loc, enrich_loc, rv_outlier_loc,
                        distribution, variant_class)
    # enrich_obj.write_rvs_w_outs_to_file(out_cut_off=2, tss_cut_off=1e4,
    #                                     af_cut_off=5e-2)
    # logger.info("Printed final set of outliers with rare variants")
    enrich_obj.loop_enrichment()
    logger.info("Completed outlier enrichment! File: " + enrich_loc)
    logger.info("All done :)")


def main():
    """Process argparse and start program."""
    parser = argparse.ArgumentParser(
        description="Associate outliers with rare variants.")
    optional = parser._action_groups.pop()
    parser.add_argument("--version", action="version", version="%(prog)s 0.1")
    # Arguments for file locations
    required = parser.add_argument_group('Required arguments')
    required.add_argument("-v", "--vcf", help="Location of VCF file",
                          required=True)
    required.add_argument("-b", "--bed", help="Gene expression file location",
                          required=True)
    optional_files = parser.add_argument_group('Optional file locations')
    optional_files.add_argument("-o", "--output", help="Output prefix")
    optional_files.add_argument("--outlier_output", help="Outlier filename")
    optional_files.add_argument("--enrich_file", help="Output file for " +
                                "enrichment odds ratios and p-values")
    # Arguments for expression outliers
    opt_out_args = parser.add_argument_group('Optional outlier arguments')
    opt_out_args.add_argument("-e", "--extrema", default=False,
                              action="store_true",
                              help="Only the most extreme value is an outlier")
    opt_out_args.add_argument("-d", "--distribution",
                              help="Outlier distribution",
                              choices=["normal", "rank", "custom"],
                              default="normal")
    opt_out_args.add_argument("-t", "--threshold", help="Expression " +
                              "threshold for defining outliers. Must be " +
                              "greater than 0 for --distribution normal or " +
                              "(0,0.5) non-inclusive with --distribution " +
                              "rank. Ignored with --distribution custom",
                              type=float, default=2.0)
    opt_out_args.add_argument("--max_outliers_per_id", help="Maximum number " +
                              "of outliers per ID", type=int, default=1000)
    # Arguments for variants
    opt_var = parser.add_argument_group('Optional variant-related arguments')
    opt_var.add_argument("--af_rare", help="AF cut-off below which a variant" +
                         "is considered rare", type=float, nargs="*",
                         default=0.01)
    opt_var.add_argument("--tss_dist", help="Variants within this distance " +
                         "of the TSS are considered", type=int, nargs="*",
                         default=1e4)
    opt_var.add_argument("--upstream", help="Only variants UPstream of TSS")
    opt_var.add_argument("--downstream", help="Only vars DOWNstream of TSS")
    opt_var.add_argument("--variant_class", help="Only variants in these " +
                         "classes will be considered", default=None,
                         choices=["intronic", "intergenic", "exonic", "UTR5",
                                  "UTR3", "splicing", "upstream", "ncRNA"])
    opt_var.add_argument("--annotations", help="Annotation file locations " +
                         "passed as a comma-separated list. Only " +
                         "variants in these annotations will be considered")
    opt_var.add_argument("--af_population", help="Space-separated " +
                         "list of populations used " +
                         "for determining the maximum observed allele " +
                         "frequency. Currently using populations with " +
                         "N>1000 on GnomAD: Non-Finnish European - NFE, " +
                         "Finnish - FIN, African/African-American - AFR, " +
                         "and global - ALL\nOther options are Latino - AMR, " +
                         "Ashkenazi Jewish - ASJ, East Asian - EAS, " +
                         "other - OTH", nargs="*",
                         default="ALL NFE FIN AFR",
                         choices=["ALL", "NFE", "FIN", "AFR",
                                  "AMR", "ASJ", "EAS", "OTH"])
    # other/utilities
    # oth_args = parser.add_argument_group('Other optional arguments')
    optional.add_argument("--processes", help="Number of CPU processes",
                          type=int, default=1)
    parser._action_groups.append(optional)
    args = parser.parse_args()
    # cprof_cmd = ('associate_outliers(args.vcf, args.bed, args.output, ' +
    #              'args.outlier_output, args.extrema, args.distribution,' +
    #              'args.variant_class, args.enrich_file)')
    # OUT_FILE = ('/sc/orga/projects/chdiTrios/Felix/dna_rna/' +
    #             'wgs_pcgc_2018_01/stats_arterial.out')
    # cProfile.run(cprof_cmd, OUT_FILE)
    # time_profile = pstats.Stats(OUT_FILE)
    # time_profile.strip_dirs().sort_stats('cumulative').print_stats(10)
    associate_outliers(args.vcf, args.bed, args.output, args.outlier_output,
                       args.extrema, args.distribution, args.variant_class,
                       args.enrich_file)


if __name__ == "__main__":
    main()
