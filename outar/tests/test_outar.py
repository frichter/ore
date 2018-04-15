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
    """Test class for Outar."""

    def test_outar(self):
        """Run test for association between RVs and outliers with mock data.

        Bad input tests:
            Abnormal VCF and BED expression inputs
            Ref genome mismatch between VCF and BED
            Not gzipped
            ANNOVAR humandb is not valid

        Calculation-based tests:
            Correct number of outliers per ID
            Correct intra-cohort AF
        """
        pass
