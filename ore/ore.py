#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Prepare variants for enrichment.

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-01-20
:Copyright: 2018, Felix Richter
:License: CC BY-SA
"""


import argparse
import re

from pkg_resources import resource_filename

from .utils import initialize_logger, checkCPUcount
from .variants import Variants
from .enrichment import Enrich
from .outliers import Outliers

"""Profiling libraries:
import cProfile
import pstats
from memory_profiler import profile
"""


# @profile
def associate_outliers(args):
    """Prepare and associate variants and outliers.

    Args:
        vcf_loc (:obj:`str`): VCF location
        bed (:obj:`str`): gene expression (phenotype) location
        output_prefix (:obj:`str`): file prefix for outputs
        outlier_postfix (:obj:`str`): file ending for outlier files
        extrema (:obj:`boolean`): T/F for using most extreme outlier
        distribution (:obj:`str`): type of outlier distribution considered
        variant_class (:obj:`str`): variant functional class to test
        enrich_loc (:obj:`str`): file location of enrichment results

    Attributes:
        outlier_obj: (:obj:`Outliers`): class instance of methods and
            states relevant to RNA expression and outlier calling

    """
    if args.output:
        output_prefix = args.output
    else:
        output_prefix = re.sub("(.*/|.vcf.gz)", "", args.vcf)
    logger = initialize_logger(log_file=output_prefix + "_ore.log",
                               logAppName="ore_status")
    checkCPUcount(args.processes)
    variants_obj = Variants(args.vcf, args.bed,
                            output_prefix=output_prefix,
                            outlier_postfix=args.outlier_output,
                            use_annovar=args.annovar,
                            annovar_dir=args.annovar_dir,
                            humandb_dir=args.humandb_dir,
                            n_processes=args.processes,
                            clean_run=args.clean_run,
                            logger=logger)
    chroms_completed = variants_obj.extract_variants_from_vcf()
    logger.info("LONG format variant abstraction done for chromosomes...\n" +
                ", ".join(chroms_completed) + "\n")
    # annovar
    if args.annovar:
        chroms_completed = variants_obj.run_annovar_wrapper()
        logger.info("Executed ANNOVAR for...\n" +
                    ", ".join(chroms_completed) + "\n")
    # find closest gene
    max_tss_dist = max(args.tss_dist)
    strand_file = resource_filename('ore', 'data/gene_strand.txt.gz')
    chroms_completed = variants_obj.label_with_closest_gene(
        upstream_only=args.upstream,
        downstream_only=args.downstream,
        max_tss_dist=max_tss_dist,
        gene_strand_data=strand_file)
    logger.info("Closest gene found for the following chromosomes..\n" +
                ", ".join(chroms_completed) + "\n")
    # other annotations
    chroms_completed = variants_obj.overlap_w_annotations_wrapper()
    logger.info("Overlaps with other bed files done for...\n" +
                ", ".join(chroms_completed) + "\n")
    # finalize variants (join variants with gene expression locations)
    chroms_completed = variants_obj.finalize_variants()
    logger.info("Final set of variants done for...\n" +
                ", ".join(chroms_completed) + "\n")
    # obtain outlier dataframe (write to file)
    outlier_obj = Outliers(pheno_loc=args.bed,
                           output_prefix=output_prefix,
                           outlier_postfix=args.outlier_output,
                           extrema=args.extrema,
                           distribution=args.distribution,
                           threshold=args.threshold)
    print("Outliers initialized...")
    outlier_obj.prepare_outliers(outlier_max=args.max_outliers_per_id,
                                 vcf_id_list=variants_obj.vcf_obj.id_list)
    logger.info("Outliers prepared")
    # output final set of outliers and calculate enrichment
    rv_outlier_loc = output_prefix + "_rv_w_outliers.txt"
    enrich_obj = Enrich(variants_obj.anno_obj.final_var_loc,
                        outlier_obj.expr_outs_loc,
                        args.enrich_file,
                        rv_outlier_loc,
                        args.distribution,
                        args.variant_class,
                        variants_obj.vcf_obj.contigs)
    enrich_obj.write_rvs_w_outs_to_file(
        out_cut_off=outlier_obj.least_extr_threshold,
        tss_cut_off=max_tss_dist,
        af_cut_off=max(args.af_rare))
    logger.info("Printed final set of outliers with rare variants")
    enrich_obj.loop_enrichment(n_processes=args.processes,
                               expr_cut_off_vec=args.threshold,
                               tss_cut_off_vec=args.tss_dist,
                               af_cut_off_vec=args.af_rare)
    logger.info("Completed outlier enrichment")
    logger.info("All done :)")


def main():
    """Process argparse and start program.

    Notes:
        None is the default value for arguments without a specified default

    """
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
    opt_out_args.add_argument("--extrema", default=False,  # "-e",
                              action="store_true",
                              help="Only the most extreme value is an outlier")
    opt_out_args.add_argument("--distribution",  # "-d",
                              help="Outlier distribution",
                              choices=["normal", "rank", "custom"],
                              default="normal")
    opt_out_args.add_argument("--threshold", help="Expression " +  # "-t",
                              "threshold for defining outliers. Must be " +
                              "greater than 0 for --distribution normal or " +
                              "(0,0.5) non-inclusive with --distribution " +
                              "rank. Ignored with --distribution custom",
                              type=float, default=2.0)
    opt_out_args.add_argument("--max_outliers_per_id", help="Maximum number " +
                              "of outliers per ID", type=int, default=None)
    # Arguments for variants
    opt_var = parser.add_argument_group('Optional variant-related arguments')
    opt_var.add_argument("--af_rare", help="AF cut-off below which a variant" +
                         "is considered rare", type=float, nargs="*",
                         default=0.01)
    opt_var.add_argument("--tss_dist", help="Variants within this distance " +
                         "of the TSS are considered", type=int, nargs="*",
                         default=1e4)
    opt_var.add_argument("--upstream", default=False, action="store_true",
                         help="Only variants UPstream of TSS")
    opt_var.add_argument("--downstream", default=False, action="store_true",
                         help="Only vars DOWNstream of TSS")
    # opt_var.add_argument("--rm_low_mapping", default=False,
    #                      action="store_true", help="Remove variants in " +
    #                      "repeats, segmental duplications, low " +
    #                      "mappability regions, and Mucin/HLA genes.")
    opt_var.add_argument("--annotations", help="Annotation file locations " +
                         "passed as a comma-separated list. Only " +
                         "variants in these annotations will be considered")
    # Variant-related arguments for ANNOVAR
    opt_annovar = parser.add_argument_group(
        'Optional arguments for using ANNOVAR')
    opt_annovar.add_argument("--annovar", default=False, action="store_true",
                             help="Use ANNOVAR to specify allele " +
                             "frequencies and functional class ()")
    # opt_annovar.add_argument("--ignore_intracohort_af", default=False,
    #                          action="store_true", help="Ignore intra-" +
    #                          "cohort AF (i.e., only " +
    #                          "use population databases for AF)")
    opt_annovar.add_argument("--variant_class", help="Only variants in " +
                             "these classes will be considered", default=None,
                             choices=["intronic", "intergenic", "exonic",
                                      "UTR5", "UTR3", "splicing", "upstream",
                                      "ncRNA"])
    opt_annovar.add_argument("--annovar_dir", help="Directory of the  " +
                             "table_annovar.pl script",
                             default=__file__ + '/annovar/')
    opt_annovar.add_argument("--humandb_dir", help="Directory of ANNOVAR " +
                             "data (refGene, ensGene, and gnomad_genome)",
                             default=__file__ + '/annovar/humandb_dir/')
    # opt_var.add_argument("--af_population", help="Space-separated " +
    #                      "list of populations used " +
    #                      "for determining the maximum observed allele " +
    #                      "frequency. Currently using populations with " +
    #                      "N>1000 on GnomAD: Non-Finnish European - NFE, " +
    #                      "Finnish - FIN, African/African-American - AFR, " +
    #                      "and global - ALL\nOther options are " +
    #                      "Latino - AMR, Ashkenazi Jewish - ASJ, " +
    #                      "East Asian - EAS, other - OTH", nargs="*",
    #                      default=["ALL", "NFE", "FIN", "AFR"],
    #                      const=["ALL", "NFE", "FIN", "AFR"],
    #                      # default="ALL NFE FIN AFR",
    #                      # const="ALL NFE FIN AFR",
    #                      # choices=["ALL NFE FIN AFR AMR ASJ EAS OTH"])
    #                      choices=["ALL", "NFE", "FIN", "AFR",
    #                               "AMR", "ASJ", "EAS", "OTH"])
    # other/utilities
    # oth_args = parser.add_argument_group('Other optional arguments')
    optional.add_argument("--processes", help="Number of CPU processes",
                          type=int, default=1)
    optional.add_argument("--clean_run", help="Delete temporary " +
                          "files from the previous run",
                          default=False, action="store_true",)
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
    associate_outliers(args)


if __name__ == "__main__":
    main()
