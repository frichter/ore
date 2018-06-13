

Running ORE
===============================

..

Flowchart
~~~~~~~~~

_static/ore_flowchart.png

Pre-processing data
~~~~~~~~~~~~~~~~~~~

Required inputs

1. Gene transcription start site (TSS) location with expression as Browser Extensible Data (BED)
2. Genotypes as Variant Call Format (VCF)


ORE uses the same inputs as FastQTL_, bgzip_ with bcftools and index with tabix_:

.. code-block:: bash

    # Prepare the VCF
    bgzip test.vcf > test.vcf.gz
    tabix -p vcf test.vcf.gz
    # Prepare the expression BED file
    bgzip test.bed > test.bed.gz
    tabix -p bed test.bed.gz


Normalize_ the VCF using :code:`bcftools norm`_. Standardize every allele into its most parsimonious and left aligned form (so that appropriate allele frequencies can be obtained from population databases)

.. code-block:: bash

    bcftools norm -f human_g1k_v37.fasta -Oz test.vcf.gz -o test.norm.vcf.gz


.. _bgzip: http://www.htslib.org/doc/bgzip.html
.. _tabix: hhttp://www.htslib.org/doc/tabix.html
.. _FastQTL: fastqtl.sourceforge.net
.. _Normalize: https://genome.sph.umich.edu/wiki/Variant_Normalization
.. _bcftools norm: http://www.htslib.org/doc/bcftools.html


Download population databases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Register with ANNOVAR_
2. You will get an email for downloading ANNOVAR. Download this and unzip/untar the file using the following commands (replacing USER_KEY with what was given in the email)

.. code-block:: bash

    wget http://www.openbioinformatics.org/annovar/download/USER_KEY/annovar.latest.tar.gz
    tar -xzvf annovar.latest.tar.gz

3. Download the databases (~20 Gb) with these commands


.. code-block:: bash

    cd annovar
    perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
    perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
    perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad_genome humandb/

4. Set the database directory as an environmental variable

.. code-block:: bash

    DB_DIR="./human_db/"


.. _ANNOVAR: http://www.openbioinformatics.org/annovar/annovar_download_form.php


Run in Python virtual environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Among many other benefits, running in a virtual environment allows one to install ORE without administrator privileges (useful for working in shared scientific computing environments).


.. code-block:: bash

    # Create the virtual environment
    virtualenv venv_ore
    # Enter the environment
    source venv_ore/bin/activate
    # make sure that pip points to the virtual environment python directory by install the latest version
    pip install --upgrade pip
    # install ORE
    pip install ore
    # leave the virtual environment
    deactivate
    # re-enter the virtual environment
    source venv_ore/bin/activate



Specify parameters
~~~~~~~~~~~~~~~~~~~

    Required arguments:
      -v VCF, --vcf VCF     Location of VCF file. Must be tabixed!
      -b BED, --bed BED     Gene expression file location. Must be tabixed!

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


Run
~~~

Run ORE using the desired parameters. Currently ORE creates many temporary files that allow for faster re-running or picking up in case of a run-time crash or error.



..  Run
    ~~~~~~~~~~~~~~~~~
    Re-run with other parameters
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Plot and interpret results
    ~~~~~~~~~~~~~~~~~~~~~~~~~





