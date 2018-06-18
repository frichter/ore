#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""VCF utility functions.

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-01-20
:Copyright: 2018, Felix Richter
:License: CC BY-SA
"""


import re


class Error(Exception):
    """Base class for exceptions in this module."""

    pass


class VCFError(Error):
    """Exception raised for errors in this module.

    Attributes:
        message -- explanation of the error

    """

    def __init__(self, message):
        """Assign error explanation to object."""
        self.message = message


def get_info_af(line_dict, alt_allele_list):
    """Obtain the allele frequency from the info field."""
    try:
        vcf_af = [i[3:] for i in line_dict["INFO"].split(";") if
                  i.startswith("AF=")]
        if len(vcf_af) != 1:
            raise VCFError("No AF in INFO field for " +
                           line_dict["#CHROM"] +
                           " " + line_dict["POS"] + " so excluding")
    except VCFError:
        vcf_af = "NA"
    finally:
        vcf_af_list = vcf_af[0].split(",")
    if len(vcf_af_list) != len(alt_allele_list):
        vcf_af_list = ["NA"] * len(alt_allele_list)
    return vcf_af_list


def replace_multiGT_with_bialleleGT(id_format_dict, alt_allele_gt, line_dict):
    """Replace multi-allelic to bi-allelic genotypes.

    Replace multi-allelic GTs (e.g., 2, 3, 4, etc) with
        bi-allelic GTs (i.e., 1) if there is a match between
        the current allele you are looping over

    TODO
        check if genotype is format other than number1/number2 or ./.
        Confirm it works with phased data (i.e., |)

    """
    # loop over proband (gt_id) and genotype from FORMAT field (gt_value)
    for ft_id, ft_field_dict in id_format_dict.items():
        try:
            gt_split_list = re.split("[|/]", ft_field_dict["GT"])
            alt_allele_gt_str = str(alt_allele_gt)
            if gt_split_list[0] == alt_allele_gt_str:
                gt_split_list[0] = "1"
            elif gt_split_list[0] == ".":
                pass
            else:
                gt_split_list[0] = "0"
            if gt_split_list[1] == alt_allele_gt_str:
                gt_split_list[1] = "1"
            elif gt_split_list[1] == ".":
                pass
            else:
                gt_split_list[1] = "0"
            gt_value_bi = "/".join(gt_split_list)
            id_format_dict[ft_id]["GT"] = gt_value_bi
        except KeyError:
            # if AD is not a field. should raise error and Exit program
            raise VCFError("No GT field in FORMAT column on line " +
                           line_dict["#CHROM"] +
                           " " + line_dict["POS"] + " " + ft_id)
        except TypeError:
            pass
            # ft_field_dict is not a dict (expected with not enough info)
        except IndexError:
            # because AD does not have all alleles listed for some reason
            raise VCFError("GT field in FORMAT column does not have ref " +
                           "and alt genotype on line " +
                           line_dict["#CHROM"] + " " +
                           line_dict["POS"] + " " + ft_id)
    return id_format_dict


def calculate_depth_and_AAR(id_format_dict, line_dict):
    """Split up the FORMAT field into a list for easier filtering."""
    for ft_id, ft_field_dict in id_format_dict.items():
        try:
            ft_field_dict["AD"] = [int(i) for i in
                                   ft_field_dict["AD"].split(",")]
            ft_field_dict["New_DP"] = (ft_field_dict["AD"][0] +
                                       ft_field_dict["AD"][1])
            # use the ref+alt depth (new_dp) to calculate alt allelic ratio
            if ft_field_dict["New_DP"] > 0:
                ft_field_dict["AAR"] = (ft_field_dict["AD"][0] /
                                        ft_field_dict["New_DP"])
            else:
                ft_field_dict["AAR"] = 0
            id_format_dict[ft_id] = ft_field_dict
        except KeyError:
            # if AD is not a field. should raise error and Exit program
            raise VCFError("No AD field in FORMAT column on line " +
                           str(line_dict["#CHROM"]) +
                           " " + str(line_dict["POS"]) + " " + ft_id)
        except ValueError:
            ft_field_dict["New_DP"] = 0
            ft_field_dict["AAR"] = 0
            pass
            # could occur because comma-separated element of
            # ft_field_dict["AD"] is not int
        except TypeError:
            pass
            # ft_field_dict not a dict (expected with not enough info)
        except IndexError:
            # because AD does not have all alleles listed for some reason
            raise VCFError("AD field in FORMAT column does not have ref" +
                           "and alt depth on line " +
                           str(line_dict["#CHROM"]) +
                           " " + str(line_dict["POS"]) + " " + ft_id)
    return id_format_dict


def convert_to_012(id_format_dict, vcf_filters, line_dict):
    """Prepare 012 genotypes from x/y genotypes.

    convert 0/0, 1/0, 0/1, and 1/1 to 0, 1, 1, and 2 respectively
        use -1 if the allele fails one of the filters (GQ < 30,
        DP < 7, AAR > 0.8 or AAR < 0.2)

    Attributes:
        gt_dict_12 (:obj:`dict`): IDs (keys) mapped to the final genotype
        gq, dp, aar

    TODO
        errors
        make sure it works with phased data
        UNIT TEST cases with missing GT/AAR/etc

    """
    gq, dp, aar = vcf_filters
    # . should be gt -1
    gt_dict_12 = {}
    for gt_id, ft_field_dict in id_format_dict.items():
        try:
            gt = ft_field_dict["GT"]
            # if GT is not a field. should raise error and Exit program
            # should be checked in `replace_multi_gt_with_biallelic_gt`
            pass_filter = True
            if gt.startswith("1/0") or gt.startswith("0/1"):
                final_gt = "1"
            elif gt.startswith("1/1"):
                final_gt = "2"
            elif gt.startswith("0/0"):
                final_gt = "0"
            else:
                final_gt = gt
            # Filter for GQ >= 30
            if int(ft_field_dict["GQ"]) < gq:
                pass_filter = False
            # Filter for depth >= 7; New_DP or DP
            if int(ft_field_dict["New_DP"]) < dp:
                pass_filter = False
            # Filter for heterozygous allelic ratio [0.2, 0.8]
            if final_gt == "1":
                if (ft_field_dict["AAR"] > aar[1]) or (ft_field_dict["AAR"]
                                                       < aar[0]):
                    pass_filter = False
            if not pass_filter:
                final_gt = "-1"
        except ValueError:
            # if ft_field_dict["DP"] == '.':
            #     final_gt = "-1"  # or maybe set to "NA"
            # else:
            print(gt_id, ft_field_dict)
            raise VCFError("Unrecognized ValueError on line " +
                           str(line_dict["#CHROM"]) +
                           " " + str(line_dict["POS"]) +
                           " " + gt_id)
        except KeyError:
            print(line_dict["FORMAT"])
            print(gt_id, ft_field_dict)
            raise VCFError("No AAR/GQ/DP field in FORMAT column on line " +
                           str(line_dict["#CHROM"]) +
                           " " + str(line_dict["POS"]) + " " + gt_id)
            # if GT is not a field. should raise error and Exit program
        except TypeError:
            # print(gt_id, ft_field_dict)
            if isinstance(ft_field_dict, str):
                final_gt = ft_field_dict
            else:
                print(gt_id, ft_field_dict)
                raise VCFError("Unrecognized TypeError on line " +
                               str(line_dict["#CHROM"]) +
                               " " + str(line_dict["POS"]) +
                               " " + gt_id)
        finally:
            gt_dict_12[gt_id] = final_gt
    return gt_dict_12
