#!/usr/bin/env nextflow

Channel
    .fromPath(params.NMF_parameters)
    .splitCsv(header:false)
    .set{ NMF_parameters }

Channel
    .fromPath(params.alterations)
    .splitCsv(header:false)
    .set{ alterations }

seed = Channel.from(1..params.iterations)


process run_random_forest {

    time = { (params.time).min }
    memory = { (params.memory + 500*(task.attempt-1)).MB }

    input:
    val input_dir_path from params.input_dir_path
    set NMF_parameters,alteration,seed from NMF_parameters.combine(alterations).combine(seed)

    output:
    path "res_*.csv" into res

    """
	#!/bin/bash

	conda activate my_base

	python $PWD/2_run_random_forest.py -i ${input_dir_path} \
                                       -p ${NMF_parameters} \
                                       -g ${alteration} \
                                       -s ${seed}
    """
}

res
    .collectFile(name: 'rf_outputs/res.tsv', keepHeader: true)
    .println { "Random forest results saved in rf_outputs/res.tsv" }
