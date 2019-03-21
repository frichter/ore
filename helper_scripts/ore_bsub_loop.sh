#BSUB -W 10:00
#BSUB -q alloc
#BUSB -n 4
#BSUB -R "rusage[mem=30000]"
#BSUB -P acc_chdiTrios
#BSUB -J atrial_intergenic
#BSUB -m mothra
#BSUB -o atrial_intergenic.stdout
#BSUB -e atrial_intergenic.stderr

module purge
module load bedtools/2.27.0 samtools/1.3 bcftools/1.6
module load python/3.5.0 py_packages/3.5
source ~/venv_ore/bin/activate

cd /sc/orga/projects/chdiTrios/Felix/dna_rna/ore/helper_scripts/
python ore_for_bsub.py

# bsub < ore_bsub_loop.sh