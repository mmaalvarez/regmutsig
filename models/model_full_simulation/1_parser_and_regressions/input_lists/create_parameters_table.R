library(tidyverse)

data.frame("number_features" = c(2,5,8,11),
           "number_tumors_per_feature" = c(5,15,20,25),
           "total_number_normal_samples" = c(50,100,150,200),
           "mean_fraction_of_nt_at_risk_mutated_in_target_bins" = format(c(5e-11,5e-8,2e-05,1e-3), scientific = FALSE),
           "mean_fraction_of_nt_at_risk_mutated_in_nontarget_bins" = format(c(1e-7,5e-8,5e-9,1e-10), scientific = FALSE)) %>% 
  expand.grid() %>%
  ## add some extra settings rows manually, e.g.:
  # add_row("number_features" = 8,
  #         "number_tumors_per_feature" = 20,
  #         "total_number_normal_samples" = 100,
  #         "mean_fraction_of_nt_at_risk_mutated_in_target_bins" = format(2e-05, scientific = FALSE),
  #         "mean_fraction_of_nt_at_risk_mutated_in_nontarget_bins" = format(5e-8, scientific = FALSE)) %>%
  arrange(number_features, number_tumors_per_feature, total_number_normal_samples, mean_fraction_of_nt_at_risk_mutated_in_target_bins, mean_fraction_of_nt_at_risk_mutated_in_nontarget_bins) %>% 
  write_csv("parameters_table.csv")
