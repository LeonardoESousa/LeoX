#!/bin/bash
#SBATCH --mem=50GB
#SBATCH --time=1-0
#SBATCH -N 1
#SBATCH -n 40
#SBATCH --partition=xeon40

module load Gaussian/16
export GAUSS_SCRDIR=/scratch/ledso/


bash $1