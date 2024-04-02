library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(valr) # for granges merging
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


### collect samples, append offset; bin by SBS; Poisson regression

# from command parameters
args = commandArgs(trailingOnly=TRUE)

sample_name = ifelse(interactive(),
                yes = "MSM0.10",
                no = gsub("\\[|\\]", "", args[1])) # after channeling in nextflow, the sample names are contained within brackets, so remove them

path_somatic_variation = ifelse(interactive(),
                                yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/data/muts_pass_",
                                no = args[2]) %>% 
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1)


# trinuc to be compared vs the other 95 in this iteration
trinuc_to_compare = ifelse(interactive(),
                           yes = "TTGG",
                           no = args[3]) %>% 
  str_replace(., "(.)(.)(.)", "\\1(\\2>\\3)")

trinuc_96 = c("A(C>A)A", "A(C>A)C", "A(C>A)G", "A(C>A)T", "A(C>G)A", "A(C>G)C", "A(C>G)G", "A(C>G)T", "A(C>T)A", "A(C>T)C", "A(C>T)G", "A(C>T)T", "A(T>A)A", "A(T>A)C", "A(T>A)G", "A(T>A)T", "A(T>C)A", "A(T>C)C", "A(T>C)G", "A(T>C)T", "A(T>G)A", "A(T>G)C", "A(T>G)G", "A(T>G)T", "C(C>A)A", "C(C>A)C", "C(C>A)G", "C(C>A)T", "C(C>G)A", "C(C>G)C", "C(C>G)G", "C(C>G)T", "C(C>T)A", "C(C>T)C", "C(C>T)G", "C(C>T)T", "C(T>A)A", "C(T>A)C", "C(T>A)G", "C(T>A)T", "C(T>C)A", "C(T>C)C", "C(T>C)G", "C(T>C)T", "C(T>G)A", "C(T>G)C", "C(T>G)G", "C(T>G)T", "G(C>A)A", "G(C>A)C", "G(C>A)G", "G(C>A)T", "G(C>G)A", "G(C>G)C", "G(C>G)G", "G(C>G)T", "G(C>T)A", "G(C>T)C", "G(C>T)G", "G(C>T)T", "G(T>A)A", "G(T>A)C", "G(T>A)G", "G(T>A)T", "G(T>C)A", "G(T>C)C", "G(T>C)G", "G(T>C)T", "G(T>G)A", "G(T>G)C", "G(T>G)G", "G(T>G)T", "T(C>A)A", "T(C>A)C", "T(C>A)G", "T(C>A)T", "T(C>G)A", "T(C>G)C", "T(C>G)G", "T(C>G)T", "T(C>T)A", "T(C>T)C", "T(C>T)G", "T(C>T)T", "T(T>A)A", "T(T>A)C", "T(T>A)G", "T(T>A)T", "T(T>C)A", "T(T>C)C", "T(T>C)G", "T(T>C)T", "T(T>G)A", "T(T>G)C", "T(T>G)G", "T(T>G)T")

# load offset of full good mappability regions (no trinuc matching, this gives us the nt_at_risk per trinuc)
offset = ifelse(interactive(),
                yes = Sys.glob("../work/*/*/good_mappability_regions_n_trinuc32_at_risk.tsv")[1],
                no = args[4]) %>% 
  read_tsv

## to keep SNVs in good mappability regions
good_mappability_regions = ifelse(interactive(),
                                  yes = "/g/strcombio/fsupek_cancer3/malvarez/chromatin_info/good_mappability_regions/k50.umap.parsed.bed",
                                  no = args[5]) %>%
  import.bed() %>% data.frame %>%
  select(seqnames, start, end) %>% 
  rename("chrom" = "seqnames")


## load sample (somatic mutations)

# find path+sample name that exists
existing_file = c()
for(file_exists in paste0(path_somatic_variation, sample_name, ".csv")){
  if(file.exists(file_exists)){
    existing_file = c(existing_file, file_exists)
  }
}

# raise error if no samples were found in any path, or if >=2 samples with same name exist in different paths (i.e. diff. datasets)
if(length(existing_file) != 1){
  stop(paste0("ERROR! No samples were found in any path OR multiple samples with the same name exist in different paths:", "\n", existing_file))
}

# load it
reg_table = read_csv(existing_file) %>%
  select(chr, start, end, tri) %>% 
  rename("chrom" = "chr") %>%
  mutate(mut_id = paste0("mut_", row_number())) %>% 
  # keep SNVs in good mappability regions
  bed_intersect(good_mappability_regions, suffix = c("_SNVs", "_k50")) %>% 
  select(c("chrom", contains("_SNVs"))) %>% 
  rename_all(~str_replace_all(., "_SNVs|_k50", "")) %>%
  select(tri) %>% 
  table %>%
  as.data.frame

## in case there are no mutations left after mapping to good_mappability_regions
if(nrow(reg_table)!=0){
  
  reg_table = reg_table %>%
    rename("mutcount" = "Freq",
           "SBS96" = ".") %>% 
    ## add offset
    merge(offset, all = T) %>%
    replace_na(list(mutcount = 0)) %>% 
    # pool all the other 95 trinuc mut types
    mutate(SBS96 = ifelse(SBS96 != trinuc_to_compare,
                          "95_other_trinucs",
                          trinuc_to_compare)) %>% 
    group_by(SBS96) %>% 
    # WARNING: the trinucAtRisk of 95_other_trinucs include the nt at risk of the trinuc_to_compare × 2 (e.g. 924386 × 2, of "ACA" for A(C>A)A) because ACA could yield also A(C>G)A and A(C>T)A; also because of this, the nt at risk for any other trinuc32 (e.g. "ATA") will be counted × 3 in the nt_at_risk, one per each trinuc96 that it can yield; basically, see the offset df, how every trinuc32 is triplicated because it yields 3 diff mutations, so every trinuc32 counts is also triplicated, and the 3 values are pooled in the final nt at risk of the 95_other_trinucs (only ×2 for the trinuc_to_compare)
    summarise(trinucAtRisk = sum(n_trinuc32_at_risk),
              trinuc_mut = sum(mutcount)) %>%
    ungroup()

  ### trinuc (vs all others) Poisson_regression
  
  # GLM for Poisson
  cat(sprintf('Running Poisson regression for %s (vs. remaining SBS95) in sample %s...\n', trinuc_to_compare, sample_name))

  y_tidy = suppressWarnings(glm(formula = "trinuc_mut ~ SBS96 + offset(log(trinucAtRisk))", 
                                family = poisson(),
                                data = reg_table)) %>% 
    ## parse output
    tidy %>% 
    # keep just the Poisson regression term
    filter(!str_detect(term, "Intercept")) %>%
    select(c(term, estimate, std.error, p.value)) %>% 
    pivot_wider(names_from = term, values_from = c(estimate, std.error, p.value)) %>%
    mutate(sample_name = sample_name) %>% 
    relocate(sample_name)
} else {
  # since there are no mutations, create an empty "y_tidy"
  
  cat(sprintf('WARNING: sample %s has 0 %s: can not run regression...\n', sample_name, trinuc_to_compare))
  
  y_tidy = data.frame(matrix(ncol = 3, nrow = 1)) %>% 
    mutate(sample_name = sample_name) %>% 
    relocate(sample_name) %>% 
    `colnames<-`(c("sample_name", paste0(c(paste0("estimate_SBS96", trinuc_to_compare), paste0("std.error_SBS96", trinuc_to_compare), paste0("p.value_SBS96", trinuc_to_compare)))))
}
gc()

## parse results
results_real_sample_SBS96_regression = y_tidy %>%
  pivot_longer(col = -matches("sample_name"), names_to = "stat_feature", values_to = "val") %>% 
  separate(stat_feature, into = c("stat", "non_ref_bin_feature_level"), sep = "_") %>% 
  mutate(feature_type = "SBS96",
         feature_name = trinuc_to_compare,
         non_ref_bin_feature_level = trinuc_to_compare) %>%
  relocate(feature_type, .after = "sample_name") %>% 
  relocate(feature_name, .after = "feature_type") %>% 
  relocate(non_ref_bin_feature_level, .after = "feature_name")
gc()

write_tsv(results_real_sample_SBS96_regression, 
          paste0("results_regression__SBS96_", trinuc_to_compare, "__sample_", sample_name, ".tsv"))
