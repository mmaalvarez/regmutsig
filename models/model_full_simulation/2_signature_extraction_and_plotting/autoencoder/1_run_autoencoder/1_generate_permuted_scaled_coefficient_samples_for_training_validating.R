library(tidyverse)
library(data.table)
library(msm) # rtnorm
library(scales)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("Position", "ggplot2")


# load coefficient tables
coefficient_table = read_tsv("../../../1_parser_and_regressions/res/results_regressions_basic_simulations.tsv") %>%
  ## extract specific parameters
  filter(total_number_features == 8 & 
           number_tumors_per_feature == 25 & 
           total_number_normal_samples == 50 & 
           mean_fraction_of_nt_at_risk_mutated_in_target_bins == 2e-05 &
           mean_fraction_of_nt_at_risk_mutated_in_nontarget_bins <= 1e-10) %>% 
  select(sample_id, factor, stat, value) %>% 
  pivot_wider(names_from = stat, values_from = value)


# write the coefficients for the final processing with the trained autoencoder (to get signatures)
original_scaled = coefficient_table %>% 
  select(!contains("conf")) %>% 
  group_by(factor) %>% 
  # scale the CI-permuted coefficients matrices to [-1, +1], because final decoding layer uses tanh as activation
  mutate(scaled_estimate = rescale(estimate,
                                   to = c(-1, 1))) %>% 
  ungroup %>% 
  select(-estimate) %>% 
  pivot_wider(names_from = factor, values_from = scaled_estimate)

# sample names must be in the first column, all other columns must be the numerical features:
#   sample_short	feature_1	feature_2	any_name ...
#     sample1	      -4	        3.4	      4
#     any_name	     4	        3.2	      0

dir.create("autoencoder_input")
write_tsv(original_scaled, "autoencoder_input/original_scaled.tsv")



##############################################################################

#### Run coefficient permutations (based on CI 95%) to generate samples for autoencoder training + validation

## Parameters and initializing of some objects

# number of permuted-samples tables (1 per epoch in autoencoder)
epochs = 10
# number of permuted samples per table
totalNumIters = 100

coefficient_Resamp = list()
set.seed(1)


## Function to generate matrices resampling betas from their CI95% distributions (instead of UPmultinomial)
resample_from_CI_and_scale = function(coefficient_table){
  
  coefficient_table %>% 
    group_by(sample_id, factor) %>% 
    summarise(resampled_estimate = rtnorm(n = 1,
                                          mean = estimate,
                                          sd = 1,
                                          lower = conf.low,
                                          upper = conf.high)) %>% 
    # scale the CI-permuted coefficients matrices to [-1, +1] (e.g. python minmax scaler), because final decoding layer uses tanh as activation
    select(!contains("conf")) %>% 
    group_by(factor) %>% 
    mutate(scaled_estimate = rescale(resampled_estimate,
                                     to = c(-1, 1))) %>% 
    ungroup %>% 
    select(-resampled_estimate) %>% 
    pivot_wider(names_from = factor, values_from = scaled_estimate)
}

# run permutations
for(epoch in 1:epochs){
  
  for(nIter in 1:totalNumIters) {
    
    print(paste0("Generating permuted sample ", nIter, " for epoch ", epoch))
    
    # for each sample (row) resample coefficients from CI95% distrs, and scale
    coefficient_Resamp[[nIter]] = resample_from_CI_and_scale(coefficient_table) %>% 
      mutate("nIter" = nIter)
    
    gc()
  }
  
  # for a given epoch, bind all permuted samples as a autoencoder training+validation data input
  coefficient_Resamp = bind_rows(coefficient_Resamp) %>% 
    # shuffle all samples
    slice(sample(1:n()))
  
  write_tsv(coefficient_Resamp, 
            paste0("autoencoder_input/permuted_coefficients_", totalNumIters, "__epoch_", epoch, ".tsv"))
  coefficient_Resamp = list()
  gc()
}

