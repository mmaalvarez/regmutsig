library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")

dir.create("rf_inputs")

# define combinations of nFact and k that we want to plot
optimal_nFact_k_list = expand.grid(seq(3,12),      # nFact
                                   seq(3,12)) %>% # K
  as.matrix %>% t %>% data.frame %>% as.list

for(optimal_nFact_k in optimal_nFact_k_list){

  nFact = optimal_nFact_k[1]
  k = optimal_nFact_k[2] 
  
  if(!identical(Sys.glob(paste0("../exposures_weights/fct", nFact, "_k", k, "_exposures.tsv")), character(0))){
    
    res_folder = paste0("rf_inputs/res_nFact-", nFact, "_k-", k)
    dir.create(res_folder)
    
    # load exposures
    exposures_table = read_tsv(paste0("../exposures_weights/fct", nFact, "_k", k, "_exposures.tsv")) %>% 
      # spread the n signature names as columns, so each Sample is in a single row with its n exposures in the same row
      pivot_wider(names_from = Signature, values_from = Exposure)
    
    #############
    ## to classify Samples between altered vs wt, for each alteration that is applied to at least n Samples
    n = 5
    alteration_at_least_n_Samples = unique(c(as.character(table(exposures_table$alteration) %>% data.frame %>% filter(Freq >= n) %>% pull(Var1))))
    
    dir.create(paste0(res_folder, "/alterations/"))
    
    for(alteration_name in alteration_at_least_n_Samples){
      
      exposures_single_alteration = exposures_table %>% 
        mutate(altered_sample = ifelse(str_detect(alteration, alteration_name), 1, 0)) %>% 
        relocate(altered_sample) %>% 
        unite("id", Sample, dataset, alteration, altered_pathway_or_treatment_type, sep = "__") %>% 
        relocate(id)
      
      write_tsv(exposures_single_alteration,
                paste0(res_folder, "/alterations/", gsub("\\ |\\-|\\/|\\,|\\.\\ ", "_", alteration_name), "_exposures.tsv"))
    }
      
    #############
    ## to classify Samples between altered vs wt, for each alteration TYPE (i.e. treatment or pathway type) that is applied to at least n Samples
    alteration_at_least_n_Samples = unique(c(as.character(table(exposures_table$altered_pathway_or_treatment_type) %>% data.frame %>% filter(Freq >= n) %>% pull(Var1))))
    
    dir.create(paste0(res_folder, "/alteration_types/"))
    
    for(alteration_type in alteration_at_least_n_Samples){
      
      exposures_single_alteration_type = exposures_table %>% 
        mutate(altered_sample = ifelse(str_detect(altered_pathway_or_treatment_type, alteration_type), 1, 0)) %>% 
        relocate(altered_sample) %>% 
        unite("id", Sample, dataset, alteration, altered_pathway_or_treatment_type, sep = "__") %>% 
        relocate(id)
      
      alteration_type = gsub("\\&", "and", alteration_type)
      
      write_tsv(exposures_single_alteration_type,
                paste0(res_folder, "/alteration_types/", gsub("\\ |\\-|\\/|\\,|\\.\\ ", "_", alteration_type), "_exposures.tsv"))
    }
  }
}
