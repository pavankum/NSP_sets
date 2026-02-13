#!/bin/bash
#SBATCH -J molprocess
#SBATCH -p standard
#SBATCH -t 4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2gb
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --mail-user=pbehara@uci.edu
#SBATCH --error=slurm_output/slurm-%j.err
#SBATCH --output=slurm_output/slurm-%j.out

source $HOME/.bashrc
conda activate fb_192


file_name=filename
cp 1.pubchem_oe_pattern_search.py ${SLURM_JOB_ID}_oe_pattern_search.py
sed -i "s/file_name/${file_name}/g" ${SLURM_JOB_ID}_oe_pattern_search.py
python ${SLURM_JOB_ID}_oe_pattern_search.py
rm ${SLURM_JOB_ID}_oe_pattern_search.py
