#!/bin/bash

conda activate nextflow

export NXF_DEFAULT_DSL=1
#export TOWER_ACCESS_TOKEN="eyJ0aWQiOiA0ODk1fS4yNDdkNzEyYmY5NGUxMmFlOTQ1OGYwYWJlYmI2MjY0YmU2Y2E4Yzdl"

mkdir -p log/

nextflow -log $PWD/log/nextflow.log run pipe.nf --features=$PWD/input_lists/features.csv \
												--sample_names=$PWD/input_lists/sample_names.csv \
												--somatic_data_paths=/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/marcel_K562/2_rm_cups_lift_to_hg19/data_hg19/muts_pass_hg19_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/data/muts_pass_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/data/muts_pass_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/petljak_2022/parse_input/res/muts_pass_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/data/muts_pass_ \
												--good_mappability_regions=/g/strcombio/fsupek_cancer3/malvarez/chromatin_info/good_mappability_regions/k50.umap.parsed.bed \
												--earliest_RT=/g/strcombio/fsupek_data/_LEGACY/epi/replicationTiming/percentileBeds2016/RepliSeq.mean8solid__eqFreqBin6of6_avgSm0_lowThr0.00.bed.gz \
												--utils=/g/strcombio/fsupek_data/users/malvarez/projects/regmutsig/bin/utils.R \
												--TMminEuclidean=0.01 \
												--TMmaxFractionRmTrinucs=0.2 \
												--flank_length=50000 \
												--time_process_a0=6 \
												--time_process_a1=2 \
												--time_process_a2=2 \
												--time_process_a3=1 \
												--time_process_b1=1 \
												--time_process_b2=0.25 \
												--memory_process_a0=22 \
												--memory_process_a1=24 \
												--memory_process_a2=7 \
												--memory_process_a3=5 \
												--memory_process_b1=5 \
												--memory_process_b2=1 \
												-resume #\	-with-tower
