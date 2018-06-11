ORE: Outlier-RV enrichment
--------------------------

Cursory use of ORE (outlier-RV enrichment) is provided here, visit the `latest ORE documentation`_ for more details. Confirm the following are installed:

    - `Python >=3.5.0`_
    - `bedtools >=2.27.0`_
    - `samtools >=1.3 and bcftools >=1.6`_

Then, on the command line, install with

.. code-block::

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


Variants and gene expression are specified with :code:`--vcf` and (line 1) :code:`--bed` (line 2), respectively. The output prefix is provided with :code:`--output` (line 3). In this example, the outlier specifications :code:`--distribution` (line 4), :code:`--threshold` (line 5), and :code:`--max_outliers_per_id` (line 6) indicate that outliers are defined using a normal distribution with a z-score more extreme than two, and samples with more than 500 outliers are excluded. Variant information is specified with :code:`--af_rare` (line 7) and :code:`--tss_dist` (line 8) to encode that variants are defined as rare with a intra-cohort allele frequency at varying thresholds (≤ 0.05, 0.01, and 0.001), and to only use variants within 5 kb of the TSS.


Usage, visit the `latest ORE documentation`_ for more

.. code-block::

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


Required arguments:
  -v VCF, --vcf VCF     Location of VCF file
  -b BED, --bed BED     Gene expression file location

Optional file locations:
  -o OUTPUT, --output OUTPUT
                        Output prefix (default is VCF prefix)
  --outlier_output OUTLIER_OUTPUT
                        Outlier filename (default is VCF prefix)
  --enrich_file ENRICH_FILE
                        Output file for enrichment odds ratios and p-values (default is VCF prefix)

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


Felix Richter <felix.richter@icahn.mssm.edu>


.. _Python >=3.5.0: https://www.python.org/downloads/release/python-350/
.. _bedtools >=2.27.0: http://bedtools.readthedocs.io/en/latest/
.. _samtools >=1.3 and bcftools >=1.6: http://www.htslib.org/download/
.. _latest ORE documentation: http://ore.readthedocs.io/en/latest/ 
