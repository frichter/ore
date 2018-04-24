.. ORE documentation master file, created by
   sphinx-quickstart on Tue Apr 24 18:53:07 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

ORE: Outlier-RV enrichment
===============================

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
~~~~~~~~~~~~~~~~~~~

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

To use ORE (outlier-RV enrichment), confirm the following are installed:

    - python >=3.5.0
    - bedtools >=2.27.0 (http://bedtools.readthedocs.io/en/latest/)
    - samtools >=1.3 and bcftools >=1.6 (both from http://www.htslib.org/download/)

Then, on the command line, install with

.. code-block:: none

    pip install ore

Example run

.. code-block:: none

    ore --vcf test.vcf.gz \
        --bed test.bed.gz \
        --enrich_file enrichment.txt \
        --distribution "normal" \
        --threshold 2 3 4 \
        --max_outliers_per_id 500 \
        --af_rare 0.05 0.01 1e-3 \
        --tss_dist 5000

Usage

.. code-block:: none

  ore [-h] [--version] -v VCF -b BED [-o OUTPUT]
           [--outlier_output OUTLIER_OUTPUT] [--enrich_file ENRICH_FILE]
           [--extrema] [--distribution {normal,rank,custom}]
           [--threshold THRESHOLD] [--max_outliers_per_id MAX_OUTLIERS_PER_ID]
           [--af_rare [AF_RARE [AF_RARE ...]]]
           [--tss_dist [TSS_DIST [TSS_DIST ...]]] [--upstream] [--downstream]
           [--annotations ANNOTATIONS] [--annovar]
           [--variant_class {intronic,intergenic,exonic,UTR5,UTR3,splicing,upstream,ncRNA}]
           [--annovar_dir ANNOVAR_DIR] [--humandb_dir HUMANDB_DIR]
           [--processes PROCESSES] [--clean_run]
