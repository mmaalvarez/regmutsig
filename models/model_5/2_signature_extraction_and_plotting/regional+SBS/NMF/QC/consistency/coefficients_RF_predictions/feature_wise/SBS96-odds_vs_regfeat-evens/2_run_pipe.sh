#!/bin/bash

# create or update the .py executable, based on the .ipynb
conda activate my_base
jupytext --set-formats ipynb,py 2_run_random_forest.ipynb --sync

conda activate nextflow

export NXF_DEFAULT_DSL=1

mkdir -p log/

nextflow -log $PWD/log/nextflow.log run $PWD/pipe.nf --input_dir_path $PWD/rf_inputs \
													 --parameters $PWD/parameters.csv \
													 --memory 1 \
													 --time 2 \
													 -resume
