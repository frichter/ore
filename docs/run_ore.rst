

Flowchart
~~~~~~~~~

_static/ore_flowchart.png

Pre-processing data
~~~~~~~~~~~~~~~~~~~

Required inputs
1. Gene transcription start site (TSS) location with expression as Browser Extensible Data (BED)
2. Genotypes as Variant Call Format (VCF)

bgzip and tabix
ORE uses the same inputs as FastQTL_, bgzip with bcftools_ and index with tabix_.
bcftools tabix -p vcf test.vcf.gz
tabix -p bed test.bed.gz

Normalize the VCF
standardize every allele into its most parsimonious and left aligned form (so that appropriate allele frequencies can be obtained from population databases)
bcftools norm test.vcf.gz

.. _bgzip: http://www.htslib.org/download/
.. _tabix: http://www.htslib.org/download/
.. _FastQTL: 


Download population databases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Register with ANNOVAR_
2. Registering will send an email with the latest ANNOVAR version. Download this and unzip/untar the file using the following commands (replacing USER_KEY with that provided in the email)

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





.. Specify parameters
    ~~~~~~~~~~~~~~~~~~
    Full list of arguments here
    Run
    ~~~~~~~~~~~~~~~~~
    Re-run with other parameters
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Plot and interpret results
    ~~~~~~~~~~~~~~~~~~~~~~~~~





