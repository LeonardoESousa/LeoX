#!/bin/bash
#SBATCH --time=2-0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=REPLACE_WITH_YOUR_PARTITION

set -e
# Replace with your cluster's Gaussian module and scratch setup.
# module load Gaussian/16
export GAUSS_SCRDIR="${TMPDIR:-/tmp}"

bash "$1"
