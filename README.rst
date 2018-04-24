ORE: Outlier-RV enrichment
--------------------------

To use ORE (outlier-RV enrichment), confirm the following are installed:

    - python >=3.5.0
    - bedtools >=2.27.0 (http://bedtools.readthedocs.io/en/latest/)
    - samtools >=1.3 and bcftools >=1.6 (both from http://www.htslib.org/download/)

Then, on the command line, install with

.. code-block::

    pip install ore

Example run

.. code-block::

    ore --vcf test.vcf.gz \
        --bed test.bed.gz \
        --enrich_file enrichment.txt \
        --distribution "normal" \
        --threshold 2 3 4 \
        --max_outliers_per_id 500 \
        --af_rare 0.05 0.01 1e-3 \
        --tss_dist 5000

Usage, visit http://ore.readthedocs.io/en/latest/ for more

.. code-block::

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

Associate outliers with rare variants.

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
                        greater than 0 for --distribution normal or (0,0.5)
                        non-inclusive with --distribution rank. Ignored with
                        --distribution custom
  --max_outliers_per_id MAX_OUTLIERS_PER_ID
                        Maximum number of outliers per ID

Optional variant-related arguments:
  --af_rare
                        AF cut-off below which a variantis considered rare
  --tss_dist
                        Variants within this distance of the TSS are
                        considered
  --upstream            Only variants UPstream of TSS
  --downstream          Only variants DOWNstream of TSS
  --annotations ANNOTATIONS
                        Annotation file locations passed as a comma-separated
                        list. Only variants in these annotations will be
                        considered

Optional arguments for using ANNOVAR:
  --annovar             Use ANNOVAR to specify allele frequencies and
                        functional class
  --variant_class
                        Only variants in these classes will be considered. Options:
                         {intronic,intergenic,exonic,UTR5,UTR3,splicing,upstream,ncRNA}
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


