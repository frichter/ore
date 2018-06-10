.. -*- coding: utf-8 -*-
.. ORE documentation master file, created by
   sphinx-quickstart on Tue Apr 24 18:53:07 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


ORE: Outlier-RV enrichment
===============================

.. toctree::
   :maxdepth: 2
   :caption: Contents:



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

    ore --vcf test.vcf.gz \
        --bed test.bed.gz \
        --enrich_file enrichment.txt \
        --distribution "normal" \
        --threshold 2 3 4 \
        --max_outliers_per_id 500 \
        --af_rare 0.05 0.01 1e-3 \
        --tss_dist 5000


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


Argument details
~~~~~~~~~~~~~~~~~~~

Required arguments:
  -v VCF, --vcf VCF     Location of VCF file
  -b BED, --bed BED     Gene expression file location

Optional file locations:
  -o OUTPUT, --output OUTPUT
                        Output prefix
  --outlier_output OUTLIER_OUTPUT
                        Outlier filename
  --enrich_file ENRICH_FILE
                        Output file for enrichment odds ratios and p-values

Optional outlier arguments:
  --extrema             Only the most extreme value is an outlier
  --distribution DISTRIBUTION
                        Outlier distribution. Options:
                        {normal,rank,custom}
  --threshold THRESHOLD
                        Expression threshold for defining outliers. Must be
                        greater than 0 for normal or (0,0.5)
                        non-inclusive with rank. Ignored with custom
  --max_outliers_per_id MAX_OUTLIERS_PER_ID
                        Maximum number of outliers per ID

Optional variant-related arguments:
  --af_rare AF_RARE
                        AF cut-off below which a variant is considered rare (space separated list e.g., 0.1 0.05)
  --af_vcf              Use the VCF AF field to define an allele as rare.
  --intracohort_rare_ac INTRACOHORT_RARE_AC
                        Allele COUNT to be used instead of intra-cohort allele
                        frequency. (still uses af_rare for population level AF
                        cut-off)
  --gq GQ
                        Minimum genotype quality each variant in each individual
  --dp DP
                        Minimum depth per variant in each individual
  --aar AAR
                        Alternate allelic ratio for heterozygous variants
                        (provide two space-separated numbers between 0 and 1,
                        e.g., 0.2 0.8)
  --tss_dist TSS_DIST
                        Variants within this distance of the TSS are
                        considered
  --upstream            Only variants UPstream of TSS
  --downstream          Only variants DOWNstream of TSS

Optional arguments for using ANNOVAR:
  --annovar             Use ANNOVAR to specify allele frequencies and
                        functional class
  --variant_class
                        Only variants in these classes will be considered. Options:
                        {intronic,intergenic,exonic,UTR5,UTR3,splicing,upstream,ncRNA}
  --exon_class
                        Only variants with these exonic impacts will be
                        considered. Options:
                        {nonsynonymous,intergenic,nonframeshift,frameshift,stopgain,stoploss}
  --refgene             Filter on RefGene function.
  --ensgene             Filter on ENSEMBL function.
  --annovar_dir ANNOVAR_DIR
                        Directory of the table_annovar.pl script
  --humandb_dir HUMANDB_DIR
                        Directory of ANNOVAR data (refGene, ensGene, and
                        gnomad_genome)

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --processes PROCESSES
                        Number of CPU processes
  --clean_run           Delete temporary files from the previous run


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