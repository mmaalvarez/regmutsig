#!/usr/bin/env nextflow



// 3 Channel systems:


// 1) in parallel each feature specified in input_lists/features.csv

Channel
    .fromPath(params.features)
    .splitCsv(header:true)
    .map{ row-> tuple(row.name, row.path) }
    .set{ features }



// 2) in parallel each sample specified in input_lists/sample_names.csv

Channel
    .fromPath(params.sample_names)
    .splitCsv(header:false)
    .set{ sample_name }

sample_name.into{sample_name_for_sample_somatic_muts_overlap_feature;sample_name_for_sample_somatic_muts_SBS96_regression}



// 3) in parallel each of the 96 trinucleotide SBS, i.e. "ACAA" == "A(C>A)A", etc.

SBS96 = Channel.from( ["ACAA", "ACAC", "ACAG", "ACAT", "ACGA", "ACGC", "ACGG", "ACGT", "ACTA", "ACTC", "ACTG", "ACTT", "ATAA", "ATAC", "ATAG", "ATAT", "ATCA", "ATCC", "ATCG", "ATCT", "ATGA", "ATGC", "ATGG", "ATGT", "CCAA", "CCAC", "CCAG", "CCAT", "CCGA", "CCGC", "CCGG", "CCGT", "CCTA", "CCTC", "CCTG", "CCTT", "CTAA", "CTAC", "CTAG", "CTAT", "CTCA", "CTCC", "CTCG", "CTCT", "CTGA", "CTGC", "CTGG", "CTGT", "GCAA", "GCAC", "GCAG", "GCAT", "GCGA", "GCGC", "GCGG", "GCGT", "GCTA", "GCTC", "GCTG", "GCTT", "GTAA", "GTAC", "GTAG", "GTAT", "GTCA", "GTCC", "GTCG", "GTCT", "GTGA", "GTGC", "GTGG", "GTGT", "TCAA", "TCAC", "TCAG", "TCAT", "TCGA", "TCGC", "TCGG", "TCGT", "TCTA", "TCTC", "TCTG", "TCTT", "TTAA", "TTAC", "TTAG", "TTAT", "TTCA", "TTCC", "TTCG", "TTCT", "TTGA", "TTGC", "TTGG", "TTGT"] )




// Run processes:


process cutoff_score {

    // calculate cutoff score of each feature, and output it so that sample_somatic_muts_overlap_feature can already start

    // takes only Channel 1

    publishDir "$PWD/res/cutoff_score_plots/", pattern: '*.score_vs_norm.mut.burdens.jpg', mode: 'copy'

    time = { (params.time_process_a0 + 5*(task.attempt-1)).hour }
    memory = { (params.memory_process_a0 + 5*(task.attempt-1)).GB }

    input:
    path utils from params.utils
    set feature_name,feature_path from features
    val good_mappability_regions from params.good_mappability_regions
    val earliest_RT from params.earliest_RT
    path sample_names from params.sample_names
    val somatic_data_paths from params.somatic_data_paths

    output:
    file('*.score_vs_norm.mut.burdens.jpg')
    tuple val(feature_name), file('*_good_mappability.tsv'), file('cutoff_score_*.tsv') into feature_file_good_mappability_and_cutoff_scores

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/a0_cutoff_score.R ${utils} ${feature_name} ${feature_path} ${good_mappability_regions} ${sample_names} ${somatic_data_paths} ${earliest_RT}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/a0_cutoff_score.R ${utils} ${feature_name} ${feature_path} ${good_mappability_regions} ${sample_names} ${somatic_data_paths} ${earliest_RT}
        fi
    else
        # no conda installed
        Rscript $PWD/bin/a0_cutoff_score.R ${utils} ${feature_name} ${feature_path} ${good_mappability_regions} ${sample_names} ${somatic_data_paths} ${earliest_RT}
    fi
    """
}



process binarize_feature_trinuc_matching_offset {

    // 1) binarize (i.e. each score becomes "high" or "low"), collapse ranges into 2 bins, add sequences and trinucleotide frequencies
    // 2) trinucleotide matching: keep the fraction of each trinuc32 that was removed in each bin, to then use to downsample the SNV counts accordingly
    // 3) calculate offset: log(# matched-trinucleotides in each genomic bin)

    // takes only Channel 1 (partially via process a0)

    time = { (params.time_process_a1 + 1*(task.attempt-1)).hour }
    memory = { (params.memory_process_a1 + 5*(task.attempt-1)).GB }

    input:
    path utils from params.utils
    set feature_name,feature_file_good_mappability,cutoff_score from feature_file_good_mappability_and_cutoff_scores
    val good_mappability_regions from params.good_mappability_regions
    val TMminEuclidean from params.TMminEuclidean
    val TMmaxFractionRmTrinucs from params.TMmaxFractionRmTrinucs
    val flank_length from params.flank_length

    output:
    tuple val(feature_name), file('feature_file_processed_*.tsv') into feature_files_processed
    tuple val(feature_name), file('trinuc_fractions_rm_per_bin_and_trinucAtRisk_*.tsv') into trinuc_fractions_rm_per_bin_and_offsets

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/a1_binarize_feature_trinuc_matching_offset.R ${utils} ${feature_name} ${feature_file_good_mappability} ${cutoff_score} ${good_mappability_regions} ${TMminEuclidean} ${TMmaxFractionRmTrinucs} ${flank_length}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/a1_binarize_feature_trinuc_matching_offset.R ${utils} ${feature_name} ${feature_file_good_mappability} ${cutoff_score} ${good_mappability_regions} ${TMminEuclidean} ${TMmaxFractionRmTrinucs} ${flank_length}
        fi
    else
        # no conda installed
        Rscript $PWD/bin/a1_binarize_feature_trinuc_matching_offset.R ${utils} ${feature_name} ${feature_file_good_mappability} ${cutoff_score} ${good_mappability_regions} ${TMminEuclidean} ${TMmaxFractionRmTrinucs} ${flank_length}
    fi
    """
}



process sample_somatic_muts_overlap_feature {

    // load sample SNVs, merge with feature map, binarize based on cutoff score for that feature map (achieved with .join)

    // combines Channel 1 (via 'feature_files_processed') and Channel 2

    time = { (params.time_process_a2 + 1*(task.attempt-1)).hour }
    memory = { (params.memory_process_a2 + 5*(task.attempt-1)).GB }

    input:
    path utils from params.utils
    set val(feature_name), file(feature_file_processed), val(sample_name) from feature_files_processed.combine(sample_name_for_sample_somatic_muts_overlap_feature)
    val somatic_data_paths from params.somatic_data_paths

    output:
    tuple val(feature_name), file('sample_*_somatic_muts_overlap_*.tsv') into samples_somatic_muts_overlap_features

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/a2_sample_somatic_muts_overlap_feature.R ${utils} ${feature_name} ${sample_name} ${somatic_data_paths} ${feature_file_processed}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/a2_sample_somatic_muts_overlap_feature.R ${utils} ${feature_name} ${sample_name} ${somatic_data_paths} ${feature_file_processed}
        fi
    else
        # no conda
        Rscript $PWD/bin/a2_sample_somatic_muts_overlap_feature.R ${utils} ${feature_name} ${sample_name} ${somatic_data_paths} ${feature_file_processed}
    fi
    """
}

// ensure that all tuples in 'samples_somatic_muts_overlap_features' that have the same 'feature_name' (1st element of the tuple), but different sample, are grouped
grouped_samples_somatic_muts_overlap_features = samples_somatic_muts_overlap_features.groupTuple(by: 0)



process apply_matching_to_somatic_muts_bind_offset_regional_feature_regression {

    // randomly remove SNVs according to matching, append offset, regional_feature regression

    // takes Channel 1 and 2, both processed by sample_somatic_muts_overlap_feature process, and .join all combinations of these channels that have a common "feature_name"

    time = { (params.time_process_a3 + 1*(task.attempt-1)).hour }
    memory = { (params.memory_process_a3 + 5*(task.attempt-1)).GB }

    input:
    set val(feature_name), file(trinuc_fractions_rm_per_bin_and_offset), file(sample_somatic_muts_overlap_feature) from trinuc_fractions_rm_per_bin_and_offsets.join(grouped_samples_somatic_muts_overlap_features)

    output:
    file 'results_regression__regional_feature_*__all_samples.tsv' into results_regional_feature_regressions

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/a3_apply_matching_to_somatic_muts_bind_offset_regional_feature_regression.R ${trinuc_fractions_rm_per_bin_and_offset} ${sample_somatic_muts_overlap_feature}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/a3_apply_matching_to_somatic_muts_bind_offset_regional_feature_regression.R ${trinuc_fractions_rm_per_bin_and_offset} ${sample_somatic_muts_overlap_feature}
        fi
    else
        # no conda
        Rscript $PWD/bin/a3_apply_matching_to_somatic_muts_bind_offset_regional_feature_regression.R ${trinuc_fractions_rm_per_bin_and_offset} ${sample_somatic_muts_overlap_feature}
    fi
    """
}

// rowbind simple regional_feature regression results of all samples - requires parsing afterwards and then manual merging to SBS regression results
results_regional_feature_regressions
    .collectFile(name: 'res/results_regional_feature_regressions_all_samples.tsv', keepHeader: true)
    .println { "Poisson regional feature regression results for all samples saved in res/results_regional_feature_regressions_all_samples.tsv" }



process offset_good_mappability_regions {

    // create offset from full good_mappability_regions, without trinuc matching, for the SBS96 regressions

    time = { (params.time_process_b1 + 1*(task.attempt-1)).hour }
    memory = { (params.memory_process_b1 + 5*(task.attempt-1)).GB }

    input:
    val good_mappability_regions from params.good_mappability_regions

    output:
    file 'good_mappability_regions_n_trinuc32_at_risk.tsv' into offset_good_mappability_regions

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/b1_offset_good_mappability_regions.R ${good_mappability_regions}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/b1_offset_good_mappability_regions.R ${good_mappability_regions}
        fi
    else
        # no conda
        Rscript $PWD/bin/b1_offset_good_mappability_regions.R ${good_mappability_regions}
    fi
    """
}



process sample_somatic_muts_SBS96_regression {

    // load sample SNVs again, combine with a SBS96, create offset, run SBS96 regression

    // takes Channel 2 and 3 and .combine them

    time = { (params.time_process_b2 + 1*(task.attempt-1)).hour }
    memory = { (params.memory_process_b2 + 5*(task.attempt-1)).GB }

    input:
    set sample_name,trinuc_to_compare from sample_name_for_sample_somatic_muts_SBS96_regression.combine(SBS96)
    val somatic_data_paths from params.somatic_data_paths
    path offset_good_mappability_regions from offset_good_mappability_regions
    val good_mappability_regions from params.good_mappability_regions

    output:
    file 'results_regression__SBS96_*__sample_*.tsv' into results_SBS96_regressions

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/b2_sample_somatic_muts_SBS96_regression.R ${sample_name} ${somatic_data_paths} ${trinuc_to_compare} ${offset_good_mappability_regions} ${good_mappability_regions} 
        else
            # no conda environment named "R"
            Rscript $PWD/bin/b2_sample_somatic_muts_SBS96_regression.R ${sample_name} ${somatic_data_paths} ${trinuc_to_compare} ${offset_good_mappability_regions} ${good_mappability_regions} 
        fi
    else
        # no conda
        Rscript $PWD/bin/b2_sample_somatic_muts_SBS96_regression.R ${sample_name} ${somatic_data_paths} ${trinuc_to_compare} ${offset_good_mappability_regions} ${good_mappability_regions} 
    fi
    """
}

// rowbind all trinuc SBS96 regression results of all samples
results_SBS96_regressions
    .collectFile(name: 'res/results_SBS96_regressions_all_samples.tsv', keepHeader: true)
    .println { "Poisson SBS96 regression results for all samples saved in res/results_SBS96_regressions_all_samples.tsv" }
