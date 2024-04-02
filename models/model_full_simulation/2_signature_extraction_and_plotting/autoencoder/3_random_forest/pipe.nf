#!/usr/bin/env nextflow

Channel
    .fromPath(params.autoencoder_parameters)
    .splitCsv(header:false)
    .set{ autoencoder_parameters }
feature = Channel.from( ['feature1','feature2','feature3','feature4','feature5','feature6','feature7','feature8'] )
test_size = Channel.from( ['0.1','0.2','0.3'] )
seed = Channel.from(1..100)

process run_random_forest {

    time = { (params.time).min }
    memory = { (params.memory + 100*(task.attempt-1)).MB }

    input:
    val input_dir_path from params.input_dir_path
    set autoencoder_parameters,feature,test_size,seed from autoencoder_parameters.combine(feature).combine(test_size).combine(seed)

    output:
    path "${autoencoder_parameters}_${feature}_${test_size}_${seed}_res.tsv" into res

    """
	#!/bin/bash

	conda activate my_base

	python $PWD/2_run_random_forest.py -i ${input_dir_path} \
                                       -p ${autoencoder_parameters} \
                                       -f ${feature} \
                                       -t ${test_size} \
                                       -s ${seed}
    """
}

res
    .collectFile(name: 'rf_outputs/res.tsv', keepHeader: true)
    .println { "Random forest results for all autoencoder_parameters, feature, test_size, and seed values concatenated and saved in rf_outputs/res.tsv" }
