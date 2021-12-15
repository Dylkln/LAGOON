#!/bin/bash

module load python/3.9
module load snakemake
module load diamond
module load conda
source /shared/projects/bluecloud/input/trophomics_pipe/conda/env/trophomics/etc/conda/activate.d/libblas_mkl_activate.sh
source activate trophomics


srun snakemake --use-conda --cores 1 --profile profile/ -j=1
