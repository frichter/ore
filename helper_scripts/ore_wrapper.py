"""Wrapper around ORE."""

import os


class OREwrapper(object):
    """Wrapper around ORE."""

    def __init__(self, home_dir, vcf, expr_f, out_class, out_prefix,
                 outlier_output, var_class, enrich_f, rm_ids):
        """Initialize the ORE wrapper."""
        self.home_dir = home_dir
        self.vcf = vcf
        self.expr_f = expr_f
        self.out_class = out_class
        self.out_prefix = out_prefix
        self.outlier_output = outlier_output
        self.var_class = var_class
        self.enrich_f = enrich_f
        self.rm_ids = rm_ids
        if out_class is 'normal':
            self.max_outs_list = [
                '100', '250', '500', '750', '1000', '2500', '5000']
        elif out_class is 'extrema':
            self.max_outs_list = [
                '100', '150', '200', '400', '600', '800', '1000']
        else:
            self.max_outs_list = [
                '100', '150', '200', '400', '600', '800', '1000']

    def run_ORE(self, sv_i, max_outs_i):
        """Run ORE."""
        if self.out_class is 'extrema':
            extrema_arg = '--extrema '
            dist_arg = 'normal'
            expr_thresh = '2'
            # max_outs_i = max_outs_list[2]
            # max_outs_arg = '--max_outliers_per_id ' + max_outs_i + ' '
            # rm_id_arg = ''
            """ Alternatively preset which IDs to remove:"""
            max_outs_i = 'custom'
            max_outs_arg = ''
            rm_id_arg = '--exclude_ids ' + self.rm_ids + ' '
        elif self.out_class is 'normal':
            extrema_arg = ''
            dist_arg = self.out_class
            expr_thresh = '2'
            # max_outs_i = max_outs_list[2]
            max_outs_arg = '--max_outliers_per_id ' + max_outs_i + ' '
            rm_id_arg = ''
            """ Alternatively preset which IDs to remove:"""
            # max_outs_i = 'custom'
            # max_outs_arg = ''
            # rm_id_arg = '--exclude_ids ' + rm_ids + ' '
        elif self.out_class is 'rank':
            extrema_arg = ''
            dist_arg = self.out_class
            expr_thresh = '0.025'
            max_outs_i = 'custom'
            max_outs_arg = ''
            rm_id_arg = '--exclude_ids ' + self.rm_ids + ' '
        expr_f_i = self.expr_f.format(sv_i)
        outlier_output_i = self.outlier_output.format(
            sv_i, self.out_class, max_outs_i)
        enrich_f_i = self.enrich_f.format(
            sv_i, self.out_class, max_outs_i, self.var_class)
        if self.var_class is 'allvars':
            var_arg = ''
        else:
            var_arg = ('--variant_class ' + self.var_class +
                       ' --refgene --ensgene ')
        ore_cmd = ('time python -m ore.ore --vcf {vcf} --bed {expr} ' +
                   '--output {out_pref} --outlier_output {outlier_pref} ' +
                   '--enrich_file {enrich} --distribution {dist} ' +
                   '--threshold {expr_thresh} ' +
                   '{extrema_arg}{max_outs_arg}{rm_id_arg}' +
                   '--af_rare 0.05 1e-2 1e-3 1e-4 1e-5 --tss_dist 5e3 1e4 ' +
                   '--annovar {var_arg}' +
                   '--humandb_dir "/sc/orga/projects/chdiTrios/whole_genome/' +
                   'humandb" --processes 3')
        ore_cmd_w_args = ore_cmd.format(
            vcf=self.vcf, expr=expr_f_i, out_pref=self.out_prefix,
            outlier_pref=outlier_output_i, enrich=enrich_f_i,
            dist=dist_arg, expr_thresh=expr_thresh,
            extrema_arg=extrema_arg,
            max_outs_arg=max_outs_arg, rm_id_arg=rm_id_arg,
            var_arg=var_arg)
        return ore_cmd_w_args

    def clean_files_after_run(self, new_data_f):
        """Move files after the run to a different directory."""
        if self.out_class is 'extrema':
            self.full_data_f = self.out_prefix + '_all_data_extrema.txt'
        elif self.out_class is 'normal':
            self.full_data_f = self.out_prefix + '_all_data.txt'
        elif self.out_class is 'rank':
            self.full_data_f = self.out_prefix + '_all_data_rank.txt'
        mv_cmd = 'mv {} {}'.format(
            self.full_data_f,
            new_data_f)
        if os.path.exists(new_data_f):
            print("Would be over-riding existing file")
        return mv_cmd

#
