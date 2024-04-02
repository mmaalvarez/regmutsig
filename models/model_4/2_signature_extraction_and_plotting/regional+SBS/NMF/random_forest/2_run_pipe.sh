#!/bin/bash


# create tables of parameters
ls rf_inputs | sed "s/res_//g" > NMF_parameters.tsv

ls rf_inputs/*/* | sort | uniq | grep exposures | sed "s/_exposures\.tsv//g" > alterations.tsv


# create or update the .py executable, based on the .ipynb
conda activate my_base
jupytext --set-formats ipynb,py 2_run_random_forest.ipynb --sync


conda activate nextflow

export NXF_DEFAULT_DSL=1

mkdir -p log/

nextflow -log $PWD/log/nextflow.log run $PWD/pipe.nf --input_dir_path $PWD/rf_inputs \
													 --NMF_parameters $PWD/NMF_parameters.tsv \
													 --alterations $PWD/alterations.tsv \
													 --iterations 50 \
													 --memory 1000 \
													 --time 2 \
													 -resume
