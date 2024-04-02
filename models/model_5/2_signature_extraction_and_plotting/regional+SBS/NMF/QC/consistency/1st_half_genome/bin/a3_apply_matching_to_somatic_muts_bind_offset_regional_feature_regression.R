library(tidyverse)
library(rlang)
library(MASS)
library(lme4)
library(broom)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("reduce", "IRanges")
conflict_prefer("expand", "tidyr")
conflict_prefer("tidy", "broom")


# from command parameters
args = commandArgs(trailingOnly=TRUE)

# load fractions of SBS96 removed in trinuc matching, and trinucAtRisk
trinuc_fractions_rm_per_bin_and_trinucAtRisk = ifelse(interactive(),
                                                      yes = Sys.glob("../work/*/*/trinuc_fractions_rm_per_bin_and_trinucAtRisk_BPDEdG.tsv")[1],
                                                      no = args[1]) %>% 
  read_tsv %>% 
  # levels as factors
  mutate_at(vars(bin),
            ~if("low" %in% unique(.)  &  "high" %in% unique(.)){
              factor(., ordered = F, levels = c('low', 'high')) # x-axis:
              # #SNVs |
              #       | ------ <-- BERdef tumor (maybe not flat, but with less negative coeff.)
              #       | \
              #       |  \  <-- BERwt tumor
              #       |___\_____ 
              #        low  high
              #        OGG1 OGG1
            }else{ ## in case levels are not "low" and "high" raise error and warn about it
              stop(paste0("The levels of at least 1 feature are not exclusively 'low' and 'high'"))
            })
feature_name = unique(trinuc_fractions_rm_per_bin_and_trinucAtRisk$feature_name)

# trinuc at risk per bin
trinucAtRisk_low_bin = trinuc_fractions_rm_per_bin_and_trinucAtRisk %>% filter(bin == "low") %>% pull(trinucAtRisk) %>% unique
trinucAtRisk_high_bin = trinuc_fractions_rm_per_bin_and_trinucAtRisk %>% filter(bin == "high") %>% pull(trinucAtRisk) %>% unique


## load sample SNVs with binarized score
sample_somatic_muts_overlap_feature = args[-1]
## if interactive
#sample_somatic_muts_overlap_feature = c(Sys.glob("../work/*/*/sample_AY9808_REP1_somatic_muts_overlap_MBD4.tsv")[1],
#                                        Sys.glob("../work/*/*/sample_MSM0.10_somatic_muts_overlap_MBD4.tsv")[1],
#                                        Sys.glob("../work/*/*/sample_MSK0.5_s4_somatic_muts_overlap_MBD4.tsv")[1])


## run each sample (if >1) independently
results_regression_all_samples = list()

for(sample_i in seq(1, length(sample_somatic_muts_overlap_feature))){
  
  sample_i_somatic_muts_overlap_feature = read_tsv(sample_somatic_muts_overlap_feature[[sample_i]])
    
  if(sum(select(sample_i_somatic_muts_overlap_feature, mutcount)) > 0){ # only parse if there are not only 0 muts
    
    sample_i_somatic_muts_overlap_feature = sample_i_somatic_muts_overlap_feature %>% 
      # levels as factors
      mutate_at(vars(bin),
                ~if(("low" %in% unique(.)  &  "high" %in% unique(.)) |
                    (unique(.) %in% c("low", "high"))){
                  factor(., ordered = F, levels = c('low', 'high'))
                }else{
                  stop(paste0("The levels of at least 1 feature are not exclusively 'low' and 'high'"))
                })
    sample_name = unique(sample_i_somatic_muts_overlap_feature$sample_name)
    
    ## create regression table
    reg_table = merge(trinuc_fractions_rm_per_bin_and_trinucAtRisk,
                      sample_i_somatic_muts_overlap_feature) %>%
      # downsample mut counts to account for the trinuc proportions removed in matching
      mutate(mutcount_downsampled = round(mutcount - mutcount * fraction_trinuc32_removed_in_bin, 0)) %>% 
      group_by(bin, trinucAtRisk) %>% 
      summarise(trinuc_mut = sum(mutcount_downsampled)) %>% 
      ungroup()
    
    ## make sure (and fix otherwise) that there are 2 bin levels
    bin_levels = unique(reg_table$bin)
    
    if(length(bin_levels) == 1){
      
      if(bin_levels == "low"){
        
        high_bin_rows = tibble(bin = factor("high", 
                                            ordered = F, levels = c("low", "high")),
                               trinucAtRisk = trinucAtRisk_high_bin,
                               trinuc_mut = 0)
        reg_table = bind_rows(reg_table,
                              high_bin_rows)
      } else {
        
        low_bin_rows = tibble(bin = factor("low", 
                                           ordered = F, levels = c("low", "high")),
                              trinucAtRisk = trinucAtRisk_low_bin,
                              trinuc_mut = 0)
        reg_table = bind_rows(reg_table,
                              low_bin_rows)
      }
    }
  } else {
    # create empty reg_table with trinuc_mut==0
    reg_table = tibble("trinuc_mut" = 0)
  }
    
  
  ### run simple Poisson regression
  
  if(sum(select(reg_table, trinuc_mut)) > 0){ # again, only do regression if there are not only 0 muts
    
    # GLM for Poisson
    cat(sprintf('Running Poisson regression for regional feature %s in sample %s...\n', feature_name, sample_name))
    y_tidy = suppressWarnings(glm(formula = "trinuc_mut ~ bin + offset(log(trinucAtRisk))", 
                                  family = poisson(),
                                  data = reg_table)) %>% 
      ## parse output
      tidy %>% 
      # keep just the pair of simple regression terms
      filter(term != "(Intercept)") %>%
      select(c(term, estimate, std.error, p.value)) %>% 
      pivot_wider(names_from = term, values_from = c(estimate, std.error, p.value)) %>%
      mutate(sample_name = sample_name) %>% 
      relocate(sample_name)
    
  } else {
    # since there are no mutations, create an empty "y_tidy"
    
    cat(sprintf('WARNING: sample %s has 0 mutations for regional feature %s: can not run regression...\n', sample_name, feature_name))
    
    y_tidy = data.frame(matrix(ncol = 3, nrow = 1)) %>% 
      mutate(sample_name = sample_name) %>% 
      relocate(sample_name) %>% 
      `colnames<-`(c("sample_name", "estimate_binhigh", "std.error_binhigh", "p.value_binhigh"))
  }
  gc()
  
  
  ## parse and add results to list
  results_regression = y_tidy %>%
    pivot_longer(col = -matches("sample_name"), names_to = "stat_feature", values_to = "val") %>% 
    separate(stat_feature, into = c("stat", "non_ref_bin_feature_level"), sep = "_") %>% 
    mutate(feature_type = "regional_feature",
           feature_name = feature_name,
           non_ref_bin_feature_level = gsub("bin", "", non_ref_bin_feature_level)) %>%
    relocate(feature_type, .after = "sample_name") %>% 
    relocate(feature_name, .after = "feature_type") %>% 
    relocate(non_ref_bin_feature_level, .after = "feature_name")
  gc()
  
  results_regression_all_samples[[sample_name]] = results_regression
}

results_regression_all_samples = bind_rows(results_regression_all_samples)

write_tsv(results_regression_all_samples, 
          paste0("results_regression__regional_feature_", feature_name, "__all_samples.tsv"))
