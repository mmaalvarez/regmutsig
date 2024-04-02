#!/bin/bash

# create or update the .py executable, based on the .ipynb
conda activate my_base
jupytext --set-formats ipynb,py autoencoder_tensorflow1.ipynb --sync

# run pipe
conda activate nextflow

export NXF_DEFAULT_DSL=1

mkdir -p log/

nextflow -log $PWD/log/nextflow.log run 2_pipe_autoencoder.nf \
	--container "/g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/containers/autoencoder/tensorflow1.15.5_gpu_jupyter_moredependencies-v1.simg" \
	--dataset_real "$PWD/autoencoder_input/original_scaled.tsv" \
	--dataset_permuted "$PWD/autoencoder_input/permuted_coefficients_*__epoch_*.tsv" \
	--memory 1 \
	-resume
