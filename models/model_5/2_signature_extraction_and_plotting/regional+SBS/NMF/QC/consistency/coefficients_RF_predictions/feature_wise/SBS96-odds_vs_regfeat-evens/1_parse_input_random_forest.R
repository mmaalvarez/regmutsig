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

regfeatures_list = read_tsv("../../../1st_half_genome/res/results_regional_feature_regressions_all_samples.tsv") %>% 
  filter(feature_type == "regional_feature") %>% 
  pull(feature_name) %>% 
  unique

critical_value = 1.28 # for CI80%; 1.96 for CI95%

CI_percs = c('CI5%', 'CI25%', 'CI50%', 'CI75%', 'CI99.99%', 'CI100%')

coefficient_tables = list()

for(genome_half in c("1st_half", "2nd_half")){
  
  coefficient_tables[[genome_half]] = list()
  
  for(feature_type in c("regional_feature", "SBS96")){
    
    coefficient_table = lapply(c(paste0("../../../", genome_half, "_genome/res/results_regional_feature_regressions_all_samples.tsv"),
                                 paste0("../../../", genome_half, "_genome/res/results_SBS96_regressions_all_samples.tsv")),
                               read_tsv) %>%
      Reduce(function(x, y) bind_rows(x, y), .) %>%
      pivot_wider(names_from = stat, values_from = val) %>% 
      # keep only either regional_features or SBS96
      filter(feature_type == .env$feature_type) %>% 
      mutate(sample_name = gsub("_REP1", "", sample_name))
    
    ### store real
    coefficient_table_original = coefficient_table %>% 
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
    
    ## store resampled with different CI%
    for(CI_perc in CI_percs){
      resampl_critical_value = case_when(CI_perc == "CI5%" ~ 0.06,
                                         CI_perc == "CI25%" ~ 0.32,
                                         CI_perc == "CI50%" ~ 0.67,
                                         CI_perc == "CI75%" ~ 1.15,
                                         CI_perc == "CI99.99%" ~ 3.89,
                                         CI_perc == "CI100%" ~ 6) # actually ~99.99999999...%
      
      coefficient_table_resampled = coefficient_table %>%    
        # for regressions that died (estimate and/or CIs >=|10|), convert the estimate and its CI to 0, as these would mean that the regression died
        mutate(conf.low = estimate - (std.error*resampl_critical_value),
               conf.high = estimate + (std.error*resampl_critical_value),
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
      
      ### if regional features, keep only 1
      if(feature_type == "regional_feature"){
        
        for(regfeature_kept in regfeatures_list){
          
          coefficient_tables[[genome_half]][[paste0(regfeature_kept, "__original")]] = coefficient_table_original %>% 
            select(id, all_of(regfeature_kept))
          
          coefficient_tables[[genome_half]][[paste0(regfeature_kept, "__resampled_", CI_perc)]] = coefficient_table_resampled %>% 
            select(id, all_of(regfeature_kept))
        }
      } else { # if SBS96 keep all 96 SBS
        coefficient_tables[[genome_half]][["SBS96__original"]] = coefficient_table_original
        coefficient_tables[[genome_half]][[paste0("SBS96__resampled_", CI_perc)]] = coefficient_table_resampled
      }
    }
  }
}


dir.create("rf_inputs")

## define combinations of predictor - topredict
all_combinations_predictor_topredict = c(map(asplit(data.frame(paste0("1st_half__", regfeatures_list, "__original"), 
                                                               paste0("2nd_half__", regfeatures_list, "__original")), 1),
                                             ~as.character(.x)),
                                         map(asplit(data.frame(paste0("1st_half__", regfeatures_list, "__original"), 
                                                               paste0("1st_half__", regfeatures_list, "__resampled_CI5%")), 1),
                                             ~as.character(.x)),
                                         map(asplit(data.frame(paste0("1st_half__", regfeatures_list, "__original"), 
                                                               paste0("1st_half__", regfeatures_list, "__resampled_CI25%")), 1),
                                             ~as.character(.x)),
                                         map(asplit(data.frame(paste0("1st_half__", regfeatures_list, "__original"), 
                                                               paste0("1st_half__", regfeatures_list, "__resampled_CI50%")), 1),
                                             ~as.character(.x)),
                                         map(asplit(data.frame(paste0("1st_half__", regfeatures_list, "__original"), 
                                                               paste0("1st_half__", regfeatures_list, "__resampled_CI75%")), 1),
                                             ~as.character(.x)),
                                         map(asplit(data.frame(paste0("1st_half__", regfeatures_list, "__original"), 
                                                               paste0("1st_half__", regfeatures_list, "__resampled_CI99.99%")), 1),
                                             ~as.character(.x)),
                                         map(asplit(data.frame(paste0("1st_half__", regfeatures_list, "__original"), 
                                                               paste0("1st_half__", regfeatures_list, "__resampled_CI100%")), 1),
                                             ~as.character(.x)),
                                         ## UPDATE: here now comparing the 1st_half__SBS96__original to the 2nd_half of regfeature, rather than to the 1st half as before
                                         map(asplit(data.frame("1st_half__SBS96__original", 
                                                               paste0("2nd_half__", regfeatures_list, "__original")), 1),
                                             ~as.character(.x)),
                                         list(c("1st_half__SBS96__original","2nd_half__SBS96__original")),
                                         list(c("1st_half__SBS96__original","1st_half__SBS96__resampled_CI5%")),
                                         list(c("1st_half__SBS96__original","1st_half__SBS96__resampled_CI25%")),
                                         list(c("1st_half__SBS96__original","1st_half__SBS96__resampled_CI50%")),
                                         list(c("1st_half__SBS96__original","1st_half__SBS96__resampled_CI75%")),
                                         list(c("1st_half__SBS96__original","1st_half__SBS96__resampled_CI99.99%")),
                                         list(c("1st_half__SBS96__original","1st_half__SBS96__resampled_CI100%")))

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


### correlation plots
coeff_list = list()
for(comparison_pair_i in seq(1, nrow(genomehalf_featuretype_predictor_and_topredict))){
  
  predictor = slice(genomehalf_featuretype_predictor_and_topredict, comparison_pair_i)$predictor
  predictor_feature_name = str_split(predictor, "__")[[1]][2]
  
  topredict = slice(genomehalf_featuretype_predictor_and_topredict, comparison_pair_i)$topredict
  topredict_feature_name = str_split(topredict, "__")[[1]][2]
  
  if(predictor_feature_name == topredict_feature_name & predictor_feature_name != "SBS96"){
    
    predictor_half = str_split(predictor, "__")[[1]][1]
    predictor_feature_type = str_split(predictor, "__")[[1]][3]
    predictor_feature_name_type = gsub(".*half__", "", predictor)
    predictor_coeff = coefficient_tables[[predictor_half]][[predictor_feature_name_type]] %>% 
      mutate(id = gsub("__.*", "", id),
             feature_name = predictor_feature_name) %>% 
      `colnames<-`(c("sample_name", "Predictor (1st genome half)", "feature_name")) %>% 
      arrange(sample_name)
    
    topredict_half = str_split(topredict, "__")[[1]][1]
    topredict_feature_type = str_split(topredict, "__")[[1]][3]
    topredict_feature_name_type = gsub(".*half__", "", topredict)
    topredict_coeff = coefficient_tables[[topredict_half]][[topredict_feature_name_type]] %>% 
      mutate(id = gsub("__.*", "", id),
             `Predicted (feature type)` = topredict_feature_type) %>% 
      `colnames<-`(c("sample_name", "Predicted", "Predicted (feature type)")) %>% 
      arrange(sample_name)
    
    coeff_list[[paste0(predictor, " vs ", topredict)]] = merge(predictor_coeff,
                                                               topredict_coeff)
  }
}
coeff_bound_tmp = bind_rows(coeff_list) %>%
  mutate(`Predicted (feature type)` = ifelse(`Predicted (feature type)` == "original",
                                             "2nd genome half", 
                                             `Predicted (feature type)`)) %>% 
  rowwise() %>% 
  mutate(`Predicted (feature type)` = ifelse(`Predicted (feature type)` != "2nd genome half",
                                             gsub("resampled_CI", "1st genome half\nresampled CI", `Predicted (feature type)`),
                                             `Predicted (feature type)`),
         `Predicted (feature type)` = gsub("CI100%", "6 SD", `Predicted (feature type)`)) %>%
  ungroup %>% 
  mutate(`Predicted (feature type)` = factor(`Predicted (feature type)`,
                                             levels = c("2nd genome half",
                                                        "1st genome half\nresampled 6 SD",
                                                        paste0("1st genome half\nresampled CI", c(99.99, 75, 50, 25, 5), "%")))) %>% 
  select(c(feature_name, `Predictor (1st genome half)`, Predicted, `Predicted (feature type)`)) %>% 
  mutate(`Predictor (1st genome half)` = ifelse(is.na(`Predictor (1st genome half)`), 0, `Predictor (1st genome half)`),
         Predicted = ifelse(is.na(Predicted), 0, Predicted))

# get Pearson per feature name and type
correlation_data = coeff_bound_tmp %>% 
  group_by(feature_name, `Predicted (feature type)`) %>% 
  summarize(Pearson = cor(`Predictor (1st genome half)`,
                           Predicted)) %>% 
  ungroup

coeff_bound = coeff_bound_tmp %>% 
  left_join(correlation_data)

pearson_corrs = ggplot(coeff_bound,
       aes(x = `Predictor (1st genome half)`,
           y = Predicted)) +
  coord_equal(ratio = 1, ylim = c(-2, 2)) +
  scale_x_continuous(n.breaks = 3,
                     labels = function(x) sprintf("%.0f", x)) +
  geom_point(size = 0.05,
             alpha = 0.5) +
  geom_text(data = coeff_bound %>% 
                      select(feature_name, `Predicted (feature type)`, Pearson) %>% 
                      distinct,
            aes(label = sprintf("r = %.1f", Pearson), 
                x = Inf, y = Inf, 
                hjust = 2.2, vjust = 1.1),
            size = 2,
            color = "blue") +
  facet_wrap(facets = vars(`Predicted (feature type)`,feature_name),
             ncol = 30,
             labeller = function (labels) {
               labels <- lapply(labels, as.character)
               list(do.call(paste, c(labels, list(sep = "\n"))))
             }) +
  theme_classic() +
  ylab("Predicted (capped to y=[-2,2])") +
  theme(text = element_text(size = 6),
        axis.text.x = element_text(size = 5), #angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 5),
        strip.text.x = element_text(size = 4), #angle = 90, hjust = 1),
        strip.background = element_blank())
ggsave("pearson_corrs.jpg",
       plot = pearson_corrs,
       device = grDevices::jpeg,
       width = 16,
       height = 9,
       dpi = 600)
