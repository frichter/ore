OutAR
-----

To use, confirm the following are installed:
module load bedtools/2.27.0
module load samtools/1.3
module load bcftools/1.6
module load python/3.5.0

## for Mac, recommend installing bedtools samtools and bcftools with
## homebrew
## confirm python3 https://docs.brew.sh/Homebrew-and-Python
brew install bedtools
brew install samtools
brew install bcftools


# alternatively (for samtools)
https://samtools.github.io/bcftools/
https://samtools.github.io/bcftools/howtos/install.html
https://www.biostars.org/p/173832/
https://gist.github.com/adefelicibus/f6fd06df1b4bb104ceeaccdd7325b856
http://www.htslib.org/
http://www.htslib.org/download/

# alternatively (for bedtools):
http://bedtools.readthedocs.io/en/latest/content/installation.html


## potential README renderer:
https://github.com/pypa/readme_renderer

# possibly tarball annovar files (depending on size)
tar -zcvf annovar_files_hg19.tar.gz annovar_files_hg19/