#!/usr/bin/env nextflow

Channel
    .fromPath(params.parameters_table)
    .splitCsv(header:true)
    .map{ row-> tuple(row.number_features, row.number_tumors_per_feature, row.total_number_normal_samples, row.mean_fraction_of_nt_at_risk_mutated_in_target_bins , row.mean_fraction_of_nt_at_risk_mutated_in_nontarget_bins) }
    .set{ parameters }

process create_simulate_samples__run_regressions {

    time = 30.min
    memory = { (params.memory + 4*(task.attempt-1)).GB }

    input:
    set number_features,number_tumors_per_feature,total_number_normal_samples,mean_fraction_of_nt_at_risk_mutated_in_target_bins,mean_fraction_of_nt_at_risk_mutated_in_nontarget_bins from parameters

    output:
    file 'results_regressions_basic_simulations.tsv' into results_regressions_basic_simulations

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/1_create_simulate_samples__run_regressions.R ${number_features} ${number_tumors_per_feature} ${total_number_normal_samples} ${mean_fraction_of_nt_at_risk_mutated_in_target_bins} ${mean_fraction_of_nt_at_risk_mutated_in_nontarget_bins}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/1_create_simulate_samples__run_regressions.R ${number_features} ${number_tumors_per_feature} ${total_number_normal_samples} ${mean_fraction_of_nt_at_risk_mutated_in_target_bins} ${mean_fraction_of_nt_at_risk_mutated_in_nontarget_bins}
        fi
    else
        # no conda
        Rscript $PWD/bin/1_create_simulate_samples__run_regressions.R ${number_features} ${number_tumors_per_feature} ${total_number_normal_samples} ${mean_fraction_of_nt_at_risk_mutated_in_target_bins} ${mean_fraction_of_nt_at_risk_mutated_in_nontarget_bins}
    fi
    """
}

// rowbind regression (simple) results of all basic simulation samples
results_regressions_basic_simulations
    .collectFile(name: 'res/results_regressions_basic_simulations.tsv', keepHeader: true)
    .println { "Simple regression results for all basic simulation samples saved in res/results_regressions_basic_simulations.tsv" }
