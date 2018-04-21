#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Test outar functions.

:Author: Felix Richter <Felix.Richter@icahn.mssm.edu>
:Date: 2018-04-15
:Copyright: 2018, Felix Richter
:License: CC BY-SA
"""

import unittest


class TestOutar(unittest.TestCase):
    """Test class for association between RVs and outliers with mock data."""

    def test_inputs(self):
        """Run tests for bad inputs.

        Abnormal VCF and BED expression inputs
        Ref genome mismatch between VCF and BED
        Not gzipped
        ANNOVAR humandb is not valid
        """
        pass

    def test_calculations(self):
        """Run tests for calculations.

        Correct number of outliers per ID
        Correct intra-cohort AF
        """
        pass

    def test_intermediate_files(self):
        """Test the intermediate files generated in a mock analysis.

        Correct number of rows and columns
        """
        pass

    def test_annovar(self):
        """Test that ANNOVAR runs correctly.

        Correct number of rows with NA allele frequencies
        """
        pass

    def test_bedtools(self):
        """Test Bedtools intersections.

        Correct number of 0s and 1s in overlap
        """
        pass

    def test_mp(self):
        """Test multiprocessing functions.

        Test that the correct number of workers were executed based on
            the number of input processes specified
        """
        pass
