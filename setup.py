#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Setuptools control."""


from setuptools import setup

# from .version import __version__

with open("README.rst", "rb") as f:
    long_descr = f.read().decode("utf-8")

# importing the version variable:
version = {}
with open("ore/version.py") as fp:
    exec(fp.read(), version)

setup(name='ore',
      packages=['ore'],
      version=version['__version__'],
      description='Associate outliers with rare variation',
      long_description=long_descr,
      entry_points={
        "console_scripts": ['ore=ore.ore:main']
        },
      url='http://github.com/frichter/ore',
      author='Felix Richter',
      author_email='felix.richter@icahn.mssm.edu',
      install_requires=[
          # 'argparse', 'logging', 'sys',
          # 'os', 'functools', 'multiprocessing', 'itertools',
          # 're', 'subprocess', 'glob', 'copy', 'gzip',
          'pandas', 'pybedtools',
          'scipy', 'matplotlib', 'numpy',
          'pysam'
      ],
      python_requires='>=3',
      # how to specify:
      # https://packaging.python.org/specifications/
      # core-metadata/#requires-external-multiple-use
      package_dir={'ore': 'ore'},
      package_data={'ore': ['data/gene_strand_hg*.txt.gz',
                            'data/hg19_genome_masks/hg19_*',
                            'annovar/annotate_variation.pl',
                            # 'annovar/humandb/*',
                            'annovar/table_annovar.pl']},
      # download_url = 'https://github.com/frichter/ore/archive/0.1.tar.gz',
      keywords=['rnaseq', 'wgs', 'outliers', 'rare_variants', 'rna', 'dna'],
      classifiers=[
        # As from http://pypi.python.org/pypi?%3Aaction=list_classifiers
        # 'Development Status :: 1 - Planning',
        # 'Development Status :: 2 - Pre-Alpha',
        'Development Status :: 3 - Alpha',
        # 'Development Status :: 4 - Beta',
        # 'Development Status :: 5 - Production/Stable',
        # 'Development Status :: 6 - Mature',
        # 'Development Status :: 7 - Inactive',
        # 'Environment :: Console',
        # 'Operating System :: MacOS :: MacOS X',
        # 'Operating System :: Microsoft :: Windows',
        # 'Operating System :: POSIX',
        # possibly specify other OS
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.'],
      # zip_safe=False,
      include_package_data=True,
      license='MIT')


"""
Sources:
- https://packaging.python.org/tutorials/distributing-packages/
- http://python-packaging.readthedocs.io/en/latest/minimal.html
- http://peterdowns.com/posts/first-time-with-pypi.html
- https://gehrcke.de/2014/02/distributing-a-python-command-line-application/
"""
