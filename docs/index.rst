.. -*- coding: utf-8 -*-
.. ORE documentation master file, created by
   sphinx-quickstart on Tue Apr 24 18:53:07 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


ORE: Outlier-RV enrichment
===============================

.. toctree::
   :maxdepth: 2
   :hidden:

   index
   run_ore




Getting started
~~~~~~~~~~~~~~~~~~~

To use ORE (outlier-RV enrichment), confirm the following are installed:

    - `Python >=3.5.0`_
    - `bedtools >=2.27.0`_
    - `samtools >=1.3 and bcftools >=1.6`_

Then, on the command line, install with pip_ using

.. code-block:: bash

    pip install ore

Example run

.. code-block:: bash
   :linenos:

    ore --vcf test.vcf.gz \
        --bed test.bed.gz \
        --output ore_results \
        --distribution normal \
        --threshold 2 3 4 \
        --max_outliers_per_id 500 \
        --af_rare 0.05 0.01 1e-3 \
        --tss_dist 5000


Variants and gene expression are specified with :code:`--vcf` (line 1) and :code:`--bed` (line 2), respectively. The output prefix is provided with :code:`--output` (line 3). In this example, the outlier specifications :code:`--distribution` (line 4), :code:`--threshold` (line 5), and :code:`--max_outliers_per_id` (line 6) indicate that outliers are defined using a normal distribution with a z-score more extreme than two, and samples with more than 500 outliers are excluded. Variant information is specified with :code:`--af_rare` (line 7) and :code:`--tss_dist` (line 8) to encode that variants are defined as rare with a intra-cohort allele frequency at varying thresholds (≤ 0.05, 0.01, and 0.001), and to only use variants within 5 kb of the TSS.



.. _Python >=3.5.0: https://www.python.org/downloads/release/python-350/
.. _bedtools >=2.27.0: http://bedtools.readthedocs.io/en/latest/
.. _samtools >=1.3 and bcftools >=1.6: http://www.htslib.org/download/
.. _pip: http://www.pip-installer.org/en/latest/index.html


Arguments
~~~~~~~~~~~~~~~~~~~


.. code-block:: none

  ore [-h] [--version] -v VCF -b BED [-o OUTPUT]
           [--outlier_output OUTLIER_OUTPUT] [--enrich_file ENRICH_FILE]
           [--extrema] [--distribution {normal,rank,custom}]
           [--threshold [THRESHOLD [THRESHOLD ...]]]
           [--max_outliers_per_id MAX_OUTLIERS_PER_ID]
           [--af_rare [AF_RARE [AF_RARE ...]]] [--af_vcf]
           [--intracohort_rare_ac INTRACOHORT_RARE_AC] [--gq GQ] [--dp DP]
           [--aar AAR AAR] [--tss_dist [TSS_DIST [TSS_DIST ...]]] [--upstream]
           [--downstream] [--annovar]
           [--variant_class {intronic,intergenic,exonic,UTR5,UTR3,splicing,upstream,ncRNA,ncRNA_exonic}]
           [--exon_class {nonsynonymous,intergenic,nonframeshift,frameshift,stopgain,stoploss}]
           [--refgene] [--ensgene] [--annovar_dir ANNOVAR_DIR]
           [--humandb_dir HUMANDB_DIR] [--processes PROCESSES] [--clean_run]


Contact
~~~~~~~~~~~~~~~~~~~

Please report **bugs** or **requests** through the `Issue
Tracker`_. Alternatively, contact Felix Richter at <felix.richter@icahn.mssm.edu>



Indices and tables
~~~~~~~~~~~~~~~~~~~

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. _Issue Tracker: https://github.com/frichter/ore/issues