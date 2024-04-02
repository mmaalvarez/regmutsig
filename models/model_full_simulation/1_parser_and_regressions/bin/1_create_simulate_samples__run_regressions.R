library(tidyverse)
library(rlang)
library(MASS)
library(broom)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("Position", "ggplot2")

args = commandArgs(trailingOnly=TRUE)


number_features = ifelse(interactive(),
                         yes = "8",
                         no = args[1]) %>% 
  as.numeric
feature_names = paste0("feature", seq(1,number_features))
features_combinations = expand.grid(rep(list(c("low", "high")), number_features))

# number tumor samples with increased mut rates PER FEATURE (highs)
number_tumors_per_feature = ifelse(interactive(),
                       yes = "20",
                       no = args[2]) %>% 
  as.numeric
tumor_sample_names = paste0("tumor", seq(1, number_tumors_per_feature))

# TOTAL number normal samples
total_number_normal_samples = ifelse(interactive(),
                        yes = 100, #number_tumors_per_feature*number_features,
                        no = args[3]) %>% 
  as.numeric
normal_sample_names = paste0("normal", seq(1, total_number_normal_samples))

# FIXED total size of the covered genome
total_nt_at_risk = 3e9
bin_size = total_nt_at_risk/nrow(features_combinations)


## mut increase per bin (sampled from a neg. binomial dist.)

# to draw random n mutations at each bin with "high" as the level of the target feature, if this is a tumor sample
# the mean is the FRACTION of the total #nucleotides that can be mutated, i.e. of the nt at risk (==bin_size)
mean_mutations_per_target_bin = ifelse(interactive(),
                      yes = "2e-05",
                      no = args[4]) %>% 
  as.numeric*bin_size

# to draw random n mutations at each bin WITHOUT "high" as the level of the target feature, if this is a tumor sample; if it is a normal sample, across all bins
# again, the mean is the FRACTION of the total #nucleotides that can be mutated, i.e. of the nt at risk (==bin_size)
# typically smaller than the ones for target bins (above)
mean_mutations_across_bins = ifelse(interactive(),
                                    yes = "5e-8",
                                    no = args[5]) %>% 
  as.numeric*bin_size



### create base regression table
base_reg_table = features_combinations %>% 
  unite(col = "bin", all_of(names(.))) %>% 
  arrange(bin) %>%
  separate(bin, into = feature_names) %>% 
  mutate(mutcount = 0,
         log_nt_at_risk = log(bin_size)) %>% 
  mutate_at(vars(contains("feature")),
            ~factor(., ordered = F, levels = c('low', 'high')))



### create simulated samples
simulated_samples = list()

for(sample_type in c("tumor", "normal")){
  
  simulated_samples[[sample_type]] = list()
  
  if(sample_type == "tumor"){
    
    for(feature in feature_names){
    
      for(sample_i in seq(1, number_tumors_per_feature)){
        
        sample_name = paste0(feature, "_tumor", sample_i)
        
        sample_reg_table = base_reg_table %>% 
          rowwise %>% 
          mutate(mutcount = ifelse(!!sym(feature) == "high",
                                   rnbinom(1, mu = mean_mutations_per_target_bin, size = mean_mutations_per_target_bin),
                                   rnbinom(1, mu = mean_mutations_across_bins, size = mean_mutations_across_bins)))
        
        simulated_samples[[sample_type]][[sample_name]] = sample_reg_table
      }    
    }
  } else if(sample_type == "normal"){
    
    for(sample_i in seq(1, total_number_normal_samples)){
      
      sample_name = paste0("normal", sample_i)
        
      sample_reg_table = base_reg_table %>% 
        rowwise %>% 
        mutate(mutcount = rnbinom(1, mu = mean_mutations_across_bins, size = mean_mutations_across_bins))
      
      simulated_samples[[sample_type]][[sample_name]] = sample_reg_table
    }
  }
}



#### regressions

simulated_samples_regressions = list()

for(sample_type in names(simulated_samples)){
  
  simulated_samples_regressions[[sample_type]] == list()
  
  for(sample_name in names(simulated_samples[[sample_type]])){
    
    # stick to a generalized linear model for the negative binomial family
    y = tryCatch({suppressWarnings(glm.nb(formula = paste0("mutcount ~ ",
                                                           paste(feature_names, collapse = " + "),
                                                           " + offset(log_nt_at_risk)"), 
                                          data = simulated_samples[[sample_type]][[sample_name]]))},
                 # track whether regression failed
                 error = function(e){"regression_failed"})
    
    # parse output (only if regression did not fail)
    if(y != "regression_failed"){
      
      y_tidy = tidy(y) %>%
        filter(term != "(Intercept)") %>%
        mutate(term = gsub("high$", "", term)) %>% 
        # calc CI 95% from std error
        mutate("conf.low" = estimate - `std.error`*1.96,
               "conf.high" = estimate + `std.error`*1.96) %>% 
        dplyr::select(c(term, estimate, contains("conf"))) %>% 
        pivot_wider(names_from = term, values_from = c(estimate, conf.low, conf.high)) %>%
        mutate(sample_id = sample_name,
               sample_type = sample_type,
               total_number_features = number_features,
               number_tumors_per_feature = number_tumors_per_feature,
               total_number_normal_samples = total_number_normal_samples,
               bin_size_bp = bin_size,
               mean_fraction_of_nt_at_risk_mutated_in_target_bins = mean_mutations_per_target_bin / bin_size,
               mean_fraction_of_nt_at_risk_mutated_in_nontarget_bins = mean_mutations_across_bins / bin_size,
               # theta value used, either because it was optimal or because it reached 10 iterations
               theta = y$theta)
      
      simulated_samples_regressions[[sample_type]][[sample_name]] = y_tidy
    }
  }
}

simulated_samples_regressions = map(simulated_samples_regressions, bind_rows) %>% 
  bind_rows() %>% 
  # since different channels will have different N features, make it long format so they can all be row-binded
  pivot_longer(cols = starts_with("estimate")|starts_with("conf"),
               names_to = "stat",
               values_to = "value") %>%
  separate(stat, into = c("stat", "factor"), sep = "_")

write_tsv(simulated_samples_regressions, "results_regressions_basic_simulations.tsv")
