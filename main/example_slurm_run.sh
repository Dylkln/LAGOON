#!/bin/bash

module load python/3.9
module load snakemake
module load diamond
module load conda
source your/path/to/your/env
source activate env_name


srun snakemake -k --use-conda --cores 1 --profile profile/ -j=1
