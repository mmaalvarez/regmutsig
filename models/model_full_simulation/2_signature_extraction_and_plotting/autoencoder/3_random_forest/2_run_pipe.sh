#!/bin/bash

# create table of parameters groups
ls rf_inputs > autoencoder_parameters.tsv


# create or update the .py executable, based on the .ipynb
conda activate my_base
jupytext --set-formats ipynb,py 2_run_random_forest.ipynb --sync


conda activate nextflow

export NXF_DEFAULT_DSL=1
#export TOWER_ACCESS_TOKEN="eyJ0aWQiOiA0ODk1fS4yNDdkNzEyYmY5NGUxMmFlOTQ1OGYwYWJlYmI2MjY0YmU2Y2E4Yzdl"

mkdir -p log/

nextflow -log $PWD/log/nextflow.log run $PWD/pipe.nf --input_dir_path $PWD/rf_inputs \
													 --autoencoder_parameters $PWD/autoencoder_parameters.tsv \
													 --memory 100 \
													 --time 2 \
													 -resume #\	-with-tower
