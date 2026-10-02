#!/bin/bash
#SBATCH --time=1-0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=REPLACE_WITH_YOUR_PARTITION

set -e
# Replace with your cluster's CREST/xTB module names.
# module load CREST xTB

# LeoX supplies the CPU count, memory, working directory and command script.
bash "$1"
