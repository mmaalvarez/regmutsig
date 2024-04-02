#!/bin/bash

conda activate nextflow

export NXF_DEFAULT_DSL=1
#export TOWER_ACCESS_TOKEN="eyJ0aWQiOiA0ODk1fS4yNDdkNzEyYmY5NGUxMmFlOTQ1OGYwYWJlYmI2MjY0YmU2Y2E4Yzdl"

mkdir -p log/

nextflow -log $PWD/log/nextflow.log run pipe.nf --parameters_table $PWD/input_lists/parameters_table.csv \
												--memory 1 \
												-resume #\-with-tower
