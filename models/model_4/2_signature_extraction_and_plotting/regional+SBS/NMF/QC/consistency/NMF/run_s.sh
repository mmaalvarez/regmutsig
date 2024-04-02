#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --mem=20G
#SBATCH -c 1

conda activate R

genome_half=$1

Rscript s.R $genome_half
