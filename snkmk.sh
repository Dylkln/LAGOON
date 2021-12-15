#!/bin/bash
#
#SBATCH --partition=long
#SBATCH --mem=150GB
#SBATCH --mail-user=klein.dylan@outlook.com
#SBATCH --mail-type=ALL

module load python/3.9
module load conda
module load snakemake
conda activate trophomics

srun snakemake -c1
