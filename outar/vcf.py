#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""VCF import, cleaning, and processing.

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-01-20
:Copyright: 2018, Felix Richter
:License: CC BY-SA
"""


import os
import re
import gzip

import pysam


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


class VCF(object):
    """Variant Call Format.

    This cleans and processes the VCF so that the contents can be used
    by the Variants class in variants.py

    Attributes:
        handle (:obj:`file_handle_type_here`): the gz file handle
        header (:obj:`str`): VCF header
        header_list (:obj:`list`): VCF #CHROM line split as a list
        id_list (:obj:`list`): list of IDs in VCF
        id_index_list (:obj:`list`): indices of IDs for each VCF line
        line (:obj:`str`): current line of VCF in utf-8 format
        line_count (:obj:`int`): line count after the header (0-based)
        prefix (:obj:`str`): prefix of VCF as well as intermediate directory
        vcf_file_loc (:obj:`str`): location of VCF

    """

    def __init__(self, vcf_file_loc, prefix=None):
        """Create a VCF object.

        This cleans and processes the VCF so that the contents can be used
        by the Variants class in variants.py

        Args:

        Attributes:
            prefix (:obj:`str`): prefix of VCF and intermediate directory
            vcf_file_loc (:obj:`str`): location of VCF

        Returns:

        Raises:
            :obj:`VCFError`: if VCF does not exist

        TODO:
            error if file is not tabix

        """
        if not os.path.exists(vcf_file_loc):
            raise VCFError("This file does not exist: '{}'"
                           .format(vcf_file_loc))
        self.vcf_file_loc = vcf_file_loc
        # if prefix:
        self.prefix = prefix
        # else:
        #     self.prefix = re.sub("(.*/|.vcf.gz)", "", vcf_file_loc)
        # get contigs
        tbx_handle = pysam.TabixFile(self.vcf_file_loc)
        self.contigs = tbx_handle.contigs
        tbx_handle.close()

    @staticmethod
    def prepare_vcf_per_chrom(current_chrom, vcf_loc, current_chrom_file_loc):
        """Convert VCF to allele-ID pairs in long format.

        Creates an instance of the encapsulating class and runs all the
            high level methods

        Args:
            current_chrom (:obj:`str`): current chromosome
            vcf_loc (:obj:`str`): location of vcf being analyzed
            current_chrom_file_loc (:obj:`str`): location of output file

        Returns:
            current_chrom (:obj:`str`): chromosome that was completed (possibly
                with prefix fail_)

        TODO:
            This is a high level function that can receive a lot of different
                errors. Figure out how to raise the same error after the
                except and
            Send excluded_ct (and other sanity checks) to a log file

        """
        vcf_obj = VCF(vcf_loc)
        vcf_obj.load_vcf()
        # remove parents: vcf_obj.remove_ids_wo_re_pattern(".*(-01|-02)$")
        if vcf_obj.ucsc_ref_genome:
            out_file = current_chrom_file_loc % current_chrom
        else:
            out_file = current_chrom_file_loc % ("chr" + current_chrom)
        if os.path.exists(out_file):
            # print("LONG format already done for chromosome", current_chrom)
            return "Not_rerun_" + current_chrom
        print("Starting chromosome", current_chrom)
        excluded_ct = vcf_obj.loop_over_vcf(current_chrom,
                                            current_chrom_file_loc)
        print("Total genotypes (var x ID) excluded for", current_chrom,
              "is", excluded_ct)
        return current_chrom

    def load_vcf(self):
        """Load VCF header and IDs.

        Attributes:
            handle (:obj:`file_handle_type_here`): the gz file handle
            line_count (:obj:`int`): line count after the header (0-based)

        """
        # self.handle = gzip.open(self.vcf_file_loc, 'rb')
        self.handle = gzip.open(self.vcf_file_loc, 'rb')
        self.parse_header()
        self.obtain_ids()
        self.determine_ref_genome()
        self.line_count = 0
        self.handle.close()

    def parse_header(self):
        """Parse VCF header.

        VCF headers have strict formats (link)

        Attributes:
            line (:obj:`str`): current line of VCF in utf-8 format
            header (:obj:`str`): VCF header

        Raises:
            :obj:`VCFError`: if header format is incorrect

        TODO:
            error handling
            Could improve/fix VCF header format here (if it's incorrect)

        """
        # convert the line to a usable format (since input in gz format)
        self.line = str(next(self.handle), 'utf-8')
        self.header = self.line
        # process vcf header, adding lines that start with ## to a header obj
        while self.line.startswith("##"):
            self.line = str(next(self.handle), 'utf-8')
            self.header += self.line

    def obtain_ids(self):
        """Extract IDs.

        #CHROM line of VCF contains all the IDs in tab-separated columns as
            well as the traditional VCF column name, so extract here.

        Attributes:
            header_list (:obj:`list`): VCF #CHROM line split as a list
            id_list (:obj:`list`): list of IDs in VCF
            id_index_list (:obj:`list`): indices of IDs for each VCF line

        TODO:
            error handling (confirm that VCF #CHROM line has correct values)

        Raises:
            :obj:`VCFError`: if no #CHROM between ## and first line of VCF
                (possibly just use default VCF column names)

        """
        if self.line.startswith("#CHROM"):
            line_list = self.line.strip().split("\t")
            self.header_list = line_list
            self.id_list = line_list[9:len(line_list)]
            self.id_index_list = [int(line_list.index(i))
                                  for i in self.id_list]

    def determine_ref_genome(self):
        """Determine which reference genome is being used.

        Attributes:
            ucsc_ref_genome (:obj:`boolean`): identifies version of ref genome
            genome_assembly (:obj:`str`): actual ref genome being used

        """
        line = str(next(self.handle), 'utf-8')
        # determine which ref genome is being used from BODY of VCF
        if line.startswith("chr"):
            self.ucsc_ref_genome = True
            # print("UCSC reference genome (ie chromosomes start with `chr`)")
            self.ref_assembly = "hg19"
        else:
            self.ucsc_ref_genome = False
            # print("GRCh reference genome (ie lines NOT starting with `chr`)")
            self.ref_assembly = "b37"
        # determine which ref genome is being used from HEADER of VCF
        assembly_match = re.search("assembly=(.*)>", self.header)
        if assembly_match:
            # print("Using genome assembly from VCF header")
            self.ref_assembly = assembly_match.group(1)
            # confirm that header and body match
            if self.ref_assembly.startswith("hg") and not self.ucsc_ref_genome:
                raise VCFError("Mismatch between genome assembly in " +
                               "body and header of VCF")
        else:
            print("No genome assembly found in VCF header, defaulting to",
                  self.ref_assembly)

    def loop_over_vcf(self, current_chrom, current_chrom_file_loc):
        """Loop over vcf.

        Loop over vcf and run functions on each line. Closes vcf once complete.
            Excludes alleles in structural variants (indicated with * as
            alt allele) or those without GATK PASS

        Attributes:
            alt_allele (:obj:`str`): current alternate allele
                from `alt_allele_list`
            gt_excluded_count (:obj:`int`): count the genotypes that were
                excluded (instead of printing to file)

        Raises:
            :obj:`VCFError`: if anything crashes out (should also close VCF)

        TODO:
            error handling
            make GATK PASS filter optional
            unit test f.close

        """
        tbx_handle = pysam.TabixFile(self.vcf_file_loc)
        self.declare_output_file_names(current_chrom_file_loc)
        self.gt_excluded_count = 0
        if self.ucsc_ref_genome:
            current_chrom = "chr" + current_chrom
        # longest chrom length is less than 3e8
        for line in tbx_handle.fetch(current_chrom, 1, 3e8):
            # str(, 'utf-8')
            self.parse_line(line)
            # * indicates alt allele is structural variant (skip these)
            if ('*' in self.line_dict["ALT"]) or \
                    (self.line_dict["FILTER"] != "PASS"):
                pass
            else:
                # process each alternate allele separately
                # important if file is multi-allelic per line
                for self.alt_allele in self.alt_allele_list:
                    self.parse_allele()
                    self.send_allele_ID_pair_to_file()
            # exit loop if current_chrom has been analyzed
            # print(self.line_count)
            # if self.line_count > 100:
            #     break
            # print the fraction of the file that's been read (if full
            # length is known)
            self.track_portion_done()
        tbx_handle.close()
        return self.gt_excluded_count

    def declare_output_file_names(self, current_chrom_file_loc):
        """One stop function for declaring all output file names (as strings).

        Attributes:
            current_chrom_file_loc (:obj:`str`): output file for variants
            error_file_loc (:obj:`str`): output file for variants with errors
            exclude_file_loc (:obj:`str`): output file for variants that were
                removed based on exclusion criteria

        """
        self.current_chrom_file_loc = current_chrom_file_loc
        self.error_file_loc = re.sub(".txt", "_unk_error.txt",
                                     current_chrom_file_loc)
        self.annovar_file_loc = re.sub("long_012", "annovar_input",
                                       current_chrom_file_loc)
        self.bed_file_loc = re.sub("long_012|txt", "bed",
                                   current_chrom_file_loc)

    def parse_line(self, line):
        """Parse each line of vcf.

        Parse each line of the vcf into lists, and the alternate
            alleles into a list

        Attributes:
            line_dict (:obj:`dict`): current line of VCF as a dictionary
                with `header_col` as keys
            alt_allele_list (:obj:`list`): alternate alleles on each line

        Raises:
            :obj:`VCFError`: if no alternate alleles are listed

        """
        line_list = line.strip().split("\t")    # str(, 'utf-8')
        self.line_dict = dict(zip(self.header_list, line_list))
        try:
            self.alt_allele_list = self.line_dict["ALT"].split(",")
            if len(self.alt_allele_list) == 0:
                raise VCFError("no ALT alleles so excluding line " +
                               self.line_count + self.line_dict["#CHROM"] +
                               " " + self.line_dict["POS"])
        except VCFError:
            self.line_dict["FILTER"] = "no_alt_alleles "

    def track_portion_done(self):
        """Calculate fraction done.

        Print the fraction of the VCF completed using `line_count`

        """
        self.line_count += 1
        if self.line_count % 1000 == 0:
            print("line of chromosome {}: {}".format(self.line_dict["#CHROM"],
                                                     self.line_count))
            # joint eqtl vcf: 43495000 or 45206995 or GQ30: 51840000
            # original WGS 350 vcf: 54630000
            # GTEx vcf: 21585872, for GQ30: 22456432,
            # GTEx GQ30 + other filters: 22450000
            # PCGC singletons: 25720000
            # obtained from bcftools norm output

    def parse_allele(self):
        """Parse one allele.

        For each alternate allele of a line, assign the corresponding
            GT and Format fields to each ID. Loop over the alternate
            alleles in loop_over_vcf. Replace all id GTs that match
            `alt_allele_gt` with 1, others with 0, keep . as .

        Attributes:
            alt_allele_gt (:obj:`int`): genotype of current
                `alt_allele` (e.g., 1, 2, 3, etc)
            id_format_dict (:obj:`dict`): dictionary mapping IDs
                (keys) to format fields (values)

        TODO:
            Confirm essential Keys are in ft_field_names (AD, GT, GQ)

        """
        # obtain the value of the allele genotype used in the FORMAT field
        self.alt_allele_gt = int(self.alt_allele_list.index(self.alt_allele))+1
        # keep subset corresponding to id_list
        self.id_format_dict = {k: self.line_dict[k] for k in
                               self.line_dict.keys() & set(self.id_list)}
        ft_field_names = self.line_dict["FORMAT"].split(":")
        for id_i, ft_field in self.id_format_dict.items():
            # ./. indicates there is not enough info to call variant
            if ft_field.startswith("./."):
                self.id_format_dict[id_i] = "no_info"
            else:
                ft_field_list = ft_field.split(":")
                if len(ft_field_list) == len(ft_field_names):
                    self.id_format_dict[id_i] = dict(zip(ft_field_names,
                                                         ft_field_list))
                else:
                    # unclear why this format field did
                    self.id_format_dict[id_i] = "unk_error, " + ft_field
        self.replace_multi_gt_with_biallelic_gt()
        self.calculate_depth_and_AAR()
        # convert genotypes to a single number (homo = 0 (from 0/0) or
        # 2 (from 1/1), het = 1 (from 1/0 or 0/1))
        self.convert_to_012()

    def replace_multi_gt_with_biallelic_gt(self):
        """Replace multi-allelic to bi-allelic genotypes.

        Replace multi-allelic GTs (e.g., 2, 3, 4, etc) with
            bi-allelic GTs (i.e., 1) if there is a match between
            the current allele you are looping over

        TODO
            check if genotype is format other than number1/number2 or ./.
            Confirm it works with phased data (i.e., |)

        """
        # loop over proband (gt_id) and genotype from FORMAT field (gt_value)
        for ft_id, ft_field_dict in self.id_format_dict.items():
            try:
                gt_split_list = ft_field_dict["GT"].split("/")
                alt_allele_gt_str = str(self.alt_allele_gt)
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
                self.id_format_dict[ft_id]["GT"] = gt_value_bi
            except KeyError:
                # if AD is not a field. should raise error and Exit program
                raise VCFError("No GT field in FORMAT column on line " +
                               self.line_count + self.line_dict["#CHROM"] +
                               " " + self.line_dict["POS"] + " " + ft_id)
            except TypeError:
                pass
                # ft_field_dict is not a dict (expected with not enough info)
            except IndexError:
                # because AD does not have all alleles listed for some reason
                raise VCFError("GT field in FORMAT column does not have ref" +
                               "and alt genotype on line " + self.line_count +
                               self.line_dict["#CHROM"] + " " +
                               self.line_dict["POS"] + " " + ft_id)
            # except:
            #     raise VCFError("Unknown error with line " +
            #                    self.line_count + self.line_dict["#CHROM"] +
            #                    " " + self.line_dict["POS"] + " " + ft_id)

    def calculate_depth_and_AAR(self):
        """Split up the FORMAT field into a list for easier filtering."""
        for ft_id, ft_field_dict in self.id_format_dict.items():
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
                self.id_format_dict[ft_id] = ft_field_dict
            except KeyError:
                # if AD is not a field. should raise error and Exit program
                raise VCFError("No AD field in FORMAT column on line " +
                               str(self.line_count) +
                               str(self.line_dict["#CHROM"]) +
                               " " + str(self.line_dict["POS"]) + " " + ft_id)
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
                               str(self.line_count) +
                               str(self.line_dict["#CHROM"]) +
                               " " + str(self.line_dict["POS"]) + " " + ft_id)

    def convert_to_012(self):
        """Prepare 012 genotypes from x/y genotypes.

        convert 0/0, 1/0, 0/1, and 1/1 to 0, 1, 1, and 2 respectively
            use -1 if the allele fails one of the filters (GQ < 30,
            DP < 7, AAR > 0.8 or AAR < 0.2)

        Attributes:
            gt_dict_12 (:obj:`dict`): IDs (keys) mapped to the final genotype

        TODO
            errors
            make sure it works with phased data
            UNIT TEST cases with missing GT/AAR/etc

        """
        # . should be gt -1
        self.gt_dict_12 = {}
        for gt_id, ft_field_dict in self.id_format_dict.items():
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
                if int(ft_field_dict["GQ"]) < 30:
                    pass_filter = False
                # Filter for depth >= 7; New_DP or DP
                if int(ft_field_dict["New_DP"]) < 7:
                    pass_filter = False
                # Filter for heterozygous allelic ratio [0.2, 0.8]
                if final_gt == "1":
                    if (ft_field_dict["AAR"] > 0.8) or (ft_field_dict["AAR"] <
                                                        0.2):
                        pass_filter = False
                if not pass_filter:
                    final_gt = "-1"
            except ValueError:
                # if ft_field_dict["DP"] == '.':
                #     final_gt = "-1"  # or maybe set to "NA"
                # else:
                print(gt_id, ft_field_dict)
                raise VCFError("Unrecognized ValueError on line " +
                               str(self.line_count) +
                               str(self.line_dict["#CHROM"]) +
                               " " + str(self.line_dict["POS"]) +
                               " " + gt_id)
            except KeyError:
                raise VCFError("No AAR/GQ/DP field in FORMAT column on line " +
                               str(self.line_count) +
                               str(self.line_dict["#CHROM"]) +
                               " " + str(self.line_dict["POS"]) + " " + gt_id)
                # if GT is not a field. should raise error and Exit program
            except TypeError:
                # print(gt_id, ft_field_dict)
                if isinstance(ft_field_dict, str):
                    final_gt = ft_field_dict
                else:
                    print(gt_id, ft_field_dict)
                    raise VCFError("Unrecognized TypeError on line " +
                                   str(self.line_count) +
                                   str(self.line_dict["#CHROM"]) +
                                   " " + str(self.line_dict["POS"]) +
                                   " " + gt_id)
            finally:
                self.gt_dict_12[gt_id] = final_gt

    def send_allele_ID_pair_to_file(self):
        """Write each allele-ID pair to a line in the file.

        Attributes:
            any_var_w_12_gt (:obj:`logical`): indicates if the alt was present
                in at least 1 allele. If True then allele info is sent to
                write_variant_position_to_anno_input_files for output

        """
        any_var_w_12_gt = False
        self.output_chrom = self.line_dict["#CHROM"]
        if not self.ucsc_ref_genome:
            self.output_chrom = "chr" + self.output_chrom
        for gt_id, gt in self.gt_dict_12.items():
            if gt in '1' or gt in '2':
                self.format_allele_ID_pair_output_line(gt_id, gt)
                # print(self.out_line)
                self.write_allele_ID_pair_to_file(self.current_chrom_file_loc)
                # out_line_ct = long_f.write(out_line + "\n")
                any_var_w_12_gt = True
            elif gt in '0' or gt in '-1' or gt in 'no_info':
                self.gt_excluded_count += 1
            else:
                self.format_allele_ID_pair_output_line(gt_id, gt)
                self.write_allele_ID_pair_to_file(self.error_file_loc)
        if any_var_w_12_gt:
            # send to ANNOVAR and bed file here
            self.write_allele_to_anno_input_files()

    def format_allele_ID_pair_output_line(self, gt_id, gt):
        """Create the output line."""
        # can uncomment to append QC parameters with
        out_pos = "%s.%s" % (self.output_chrom,
                             self.line_dict["POS"])
        self.out_line = "\t".join([out_pos, gt_id, gt, self.line_dict["REF"],
                                   self.alt_allele])

    def write_allele_ID_pair_to_file(self, chrom_out_file):
        """Write to output file."""
        chrom_file_current = chrom_out_file % self.output_chrom
        with open(chrom_file_current, 'a') as chrom_f:
            chrom_f.write(self.out_line + "\n")  # possibly use _ =
            # https://stackoverflow.com/questions/41149781/how-to-prevent-
            # f-write-to-output-the-number-of-characters-written

    def write_allele_to_anno_input_files(self):
        """Prepare variants for annotations.

        ANNOVAR input is 1-based: chrom | start | end | ref | alt | var_id
        Bed overlap input is 0-based: chrom | start | end | ref | alt

        """
        # annovar 1-based outputs
        var_id = ".".join([self.output_chrom,
                           self.line_dict["POS"],
                           self.line_dict["REF"],
                           self.alt_allele])

        deletion_end = str(int(self.line_dict["POS"]) +
                           len(self.line_dict["REF"]) - 1)
        annovar_out_line = "\t".join([self.output_chrom,
                                      self.line_dict["POS"],
                                      deletion_end,
                                      self.line_dict["REF"],
                                      self.alt_allele,
                                      var_id])
        annovar_out_loc = self.annovar_file_loc % self.output_chrom
        with open(annovar_out_loc, 'a') as chrom_f:
            chrom_f.write(annovar_out_line + "\n")  # possibly use _ =
        # bed file 0-based output
        zero_based_start = str(int(self.line_dict["POS"]) - 1)
        bed_out_line = "\t".join([self.output_chrom,
                                  zero_based_start,
                                  deletion_end,
                                  self.line_dict["REF"],
                                  self.alt_allele])
        bed_out_loc = self.bed_file_loc % self.output_chrom
        with open(bed_out_loc, 'a') as chrom_f:
            chrom_f.write(bed_out_line + "\n")  # possibly use _ =
