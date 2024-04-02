library(tidyverse)
library(msm) # rtnorm
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("Position", "ggplot2")


### as positive control resample 1 time from each matrix
resample_from_CI = function(coefficient_table){
  summarise(coefficient_table, resampled_estimate = rtnorm(n = 1,
                                                           mean = estimate,
                                                           sd = 1, 
                                                           lower = conf.low,
                                                           upper = conf.high))
}

set.seed(1)

coefficient_tables = list()
for(genome_half in c("1st_half", "2nd_half")){
  
  coefficient_tables[[genome_half]] = list()
  
  for(feature_type in c("regional_feature", "SBS96")){
  
    coefficient_table = lapply(c(paste0("../../", genome_half, "_genome/res/results_regional_feature_regressions_all_samples.tsv"),
                                 paste0("../../", genome_half, "_genome/res/results_SBS96_regressions_all_samples.tsv")),
                               read_tsv) %>%
      Reduce(function(x, y) bind_rows(x, y), .) %>%
      pivot_wider(names_from = stat, values_from = val) %>% 
      # keep only either regional_features or SBS96
      filter(feature_type == .env$feature_type) %>% 
      mutate(sample_name = gsub("_REP1", "", sample_name))
    
    ## store real
    critical_value = 1.28 # for CI80%; 1.96 for CI95%
    coefficient_tables[[genome_half]][[paste0(feature_type, "__original")]] = coefficient_table %>% 
      # for regressions that died (estimate and/or CIs >=|10|), convert the estimate and its CI to 0, as these would mean that the regression died
      mutate(conf.low = estimate - (std.error*critical_value),
             conf.high = estimate + (std.error*critical_value),
             reg_died = ifelse(abs(estimate) >= 10 | conf.low <= -10 | conf.high >= 10, "yes", "no"),
             across(starts_with("estimate") | starts_with("conf."),
                    ~ifelse(abs(estimate) >= 10 | conf.low <= -10 | conf.high >= 10, 0, .)),
             estimate = ifelse(is.na(estimate), 0, estimate),
             conf.low = ifelse(is.na(conf.low), 0, conf.low),
             conf.high = ifelse(is.na(conf.high), 0, conf.high)) %>%
      select(sample_name, feature_type, feature_name, estimate) %>% 
      arrange(sample_name, feature_name) %>% 
      pivot_wider(names_from = feature_name, values_from = estimate) %>% 
      mutate("genome" = genome_half) %>% 
      relocate(genome, .after = "sample_name") %>% 
      unite("id", sample_name, genome, feature_type, sep = "__")
    
    ## store resampled with different CI resamplings
    CI_percs = c('CI5%', 'CI25%', 'CI50%', 'CI75%', 'CI99.99%', 'CI100%')

    for(CI_perc in CI_percs){
    
      critical_value = case_when(CI_perc == "CI5%" ~ 0.06,
                                 CI_perc == "CI25%" ~ 0.32,
                                 CI_perc == "CI50%" ~ 0.67,
                                 CI_perc == "CI75%" ~ 1.15,
                                 CI_perc == "CI99.99%" ~ 3.89,
                                 CI_perc == "CI100%" ~ 6)
      
      coefficient_tables[[genome_half]][[paste0(feature_type, "__resampled_", CI_perc)]] = coefficient_table %>%    
        # for regressions that died (estimate and/or CIs >=|10|), convert the estimate and its CI to 0, as these would mean that the regression died
        mutate(conf.low = estimate - (std.error*critical_value),
               conf.high = estimate + (std.error*critical_value),
               reg_died = ifelse(abs(estimate) >= 10 | conf.low <= -10 | conf.high >= 10, "yes", "no"),
               across(starts_with("estimate") | starts_with("conf."),
                      ~ifelse(abs(estimate) >= 10 | conf.low <= -10 | conf.high >= 10, 0, .)),
               estimate = ifelse(is.na(estimate), 0, estimate),
               conf.low = ifelse(is.na(conf.low), 0, conf.low),
               conf.high = ifelse(is.na(conf.high), 0, conf.high)) %>%
        select(sample_name, feature_type, feature_name, estimate, conf.low, conf.high) %>% 
        arrange(sample_name, feature_name) %>% 
        mutate("genome" = genome_half) %>% 
        relocate(genome, .after = "sample_name") %>% 
        unite("id", sample_name, genome, feature_type, sep = "__") %>% 
        group_by(id, feature_name) %>% 
        resample_from_CI(.) %>% 
        ungroup %>% 
        pivot_wider(names_from = feature_name, values_from = resampled_estimate) 
      gc()
    }
  }
}


dir.create("rf_inputs")

## define combinations of predictor - topredict
# all_combinations_predictor_topredict = expand.grid(#c("1st_half","2nd_half"),
#                                                    c("1st_half__regional_feature", "2nd_half__SBS96"),
#                                                    c("original", "resampled")) %>% 
#   data.frame %>% unite("id", Var1, Var2, sep = "__") %>% 
#   pull(id) %>% 
#   combn(., 2, simplify = F)
all_combinations_predictor_topredict = list(c("1st_half__regional_feature__original","2nd_half__regional_feature__original"),
                                            c("1st_half__regional_feature__original","1st_half__regional_feature__resampled_CI5%"),
                                            c("1st_half__regional_feature__original","1st_half__regional_feature__resampled_CI25%"),
                                            c("1st_half__regional_feature__original","1st_half__regional_feature__resampled_CI50%"),
                                            c("1st_half__regional_feature__original","1st_half__regional_feature__resampled_CI75%"),
                                            c("1st_half__regional_feature__original","1st_half__regional_feature__resampled_CI99.99%"),
                                            c("1st_half__regional_feature__original","1st_half__regional_feature__resampled_CI100%"),
                                            c("1st_half__SBS96__original","1st_half__regional_feature__original"),
                                            c("1st_half__SBS96__original","2nd_half__SBS96__original"),
                                            c("1st_half__SBS96__original","1st_half__SBS96__resampled_CI5%"),
                                            c("1st_half__SBS96__original","1st_half__SBS96__resampled_CI25%"),
                                            c("1st_half__SBS96__original","1st_half__SBS96__resampled_CI50%"),
                                            c("1st_half__SBS96__original","1st_half__SBS96__resampled_CI75%"),
                                            c("1st_half__SBS96__original","1st_half__SBS96__resampled_CI99.99%"),
                                            c("1st_half__SBS96__original","1st_half__SBS96__resampled_CI100%"))

genomehalf_featuretype_predictor_and_topredict = tibble("predictor" = as.character(), 
                                                        "topredict" = as.character())
res = list()
for(predictor_topredict_pair in all_combinations_predictor_topredict){
  
  for(comparison_order in c(1)){#,2)){ # (if there is also the ',2', every table first acts as predictor, and then it acts as topredict)
    if(comparison_order == 1){
      predictor_data = predictor_topredict_pair[1]
      topredict_data = predictor_topredict_pair[2] 
    } else if(comparison_order == 2){
      predictor_data = predictor_topredict_pair[2]
      topredict_data = predictor_topredict_pair[1] 
    }
    ### parse them accordingly as predictor and topredict
    
    res_folder = paste0("rf_inputs/predictor-", predictor_data, "___topredict-", topredict_data)
    dir.create(res_folder)
    
    predictor_name = strsplit(predictor_data, split = "__")[[1]]
    topredict_name = strsplit(topredict_data, split = "__")[[1]]
    
    predictor_coeff = coefficient_tables[[predictor_name[1]]][[paste0(predictor_name[2], "__", predictor_name[3])]]
    topredict_coeff = coefficient_tables[[topredict_name[1]]][[paste0(topredict_name[2], "__", topredict_name[3])]]
      
    write_tsv(predictor_coeff,
              paste0(res_folder, "/predictor.tsv"))
    write_tsv(topredict_coeff,
              paste0(res_folder, "/topredict.tsv"))
    
    genomehalf_featuretype_predictor_and_topredict = bind_rows(genomehalf_featuretype_predictor_and_topredict,
                                                               tibble("predictor" = predictor_data, 
                                                                      "topredict" = topredict_data))
  }
}

write_delim(genomehalf_featuretype_predictor_and_topredict, "parameters.csv", delim = ",")
