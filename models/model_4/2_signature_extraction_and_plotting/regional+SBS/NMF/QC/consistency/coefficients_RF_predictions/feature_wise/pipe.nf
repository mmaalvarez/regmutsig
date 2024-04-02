#!/usr/bin/env nextflow

Channel
    .fromPath(params.parameters)
    .splitCsv(header:true)
    .map{ row-> tuple(row.predictor, row.topredict) }
    .set{ parameters }

process run_random_forest {

    time = { (params.time).min }
    memory = { (params.memory + 1*(task.attempt-1)).GB }

    input:
    val input_dir_path from params.input_dir_path
    set predictor,topredict from parameters

    output:
    path "res_*.tsv" into res

    """
	#!/bin/bash

	conda activate my_base

	python $PWD/2_run_random_forest.py -i ${input_dir_path} \
                                       --predictor ${predictor} \
                                       --topredict ${topredict}
    """
}

res
    .collectFile(name: 'rf_outputs/res.tsv', keepHeader: true)
    .println { "Random forest results saved in rf_outputs/res.tsv" }
