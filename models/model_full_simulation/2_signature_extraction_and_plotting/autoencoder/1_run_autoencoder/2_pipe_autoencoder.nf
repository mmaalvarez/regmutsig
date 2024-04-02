#!/usr/bin/env nextflow

validation = Channel.from( ['0.1','0.2'] )
num_components = Channel.from( ['7','8','9'] )
depth = Channel.from( ['1','2'] )
hidden_dim = Channel.from( ['1','2'] )
batch_size = Channel.from( ['100','200'] )
epochs = Channel.from( ['200','300'] )
activation_function = Channel.from( ['relu','tanh'] )
weight_initializer = Channel.from( ['glorot_uniform','glorot_normal'] )
loss_function = Channel.from( ['mean_squared_error','cosine_similarity'] )
learning_rate = Channel.from( ['0.0005','0.001'] )
optimizer = Channel.from( ['Adam','SGD'] )

process run_autoencoder {

    publishDir "$PWD/autoencoder_output/", mode: 'copy'

    time = 4.min
    memory = { (params.memory + 1*(task.attempt-1)).GB }

    input:
    // fixed paths
    path container from params.container
    path dataset_real from params.dataset_real
    val dataset_permuted from params.dataset_permuted
    // variables (channels)
    set validation,num_components,depth,hidden_dim,batch_size,epochs,activation_function,weight_initializer,loss_function,learning_rate,optimizer from validation.combine(num_components).combine(depth).combine(hidden_dim).combine(batch_size).combine(epochs).combine(activation_function).combine(weight_initializer).combine(loss_function).combine(learning_rate).combine(optimizer)

    output:
    path '*-validation_*-components_*-depth_*-hiddendim_*-batchsize_*-epochs_*_*_*_*-learningrate_*_*-meansumactivity_*-seed/*'

    """
    #!/usr/bin/env bash

    #conda activate singularity
    conda activate my_base

    #singularity exec ${container} python $PWD/autoencoder_tensorflow1.py
    python $PWD/autoencoder_tensorflow1.py \
        --dataset_real ${dataset_real} \
        --dataset_permuted "${dataset_permuted}" \
        --validation ${validation} \
        --num_components ${num_components} \
        --depth ${depth} \
        --hidden_dim ${hidden_dim} \
        --batch_size ${batch_size} \
        --epochs ${epochs} \
        --activation_function ${activation_function} \
        --weight_initializer ${weight_initializer} \
        --loss_function ${loss_function} \
        --learning_rate ${learning_rate} \
        --optimizer ${optimizer}
    """
}
