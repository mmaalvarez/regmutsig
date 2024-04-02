#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=2G

conda activate singularity

singularity pull docker://mvpandapaw/tensorflow1.15.5_gpu_jupyter_moredependencies:v1
