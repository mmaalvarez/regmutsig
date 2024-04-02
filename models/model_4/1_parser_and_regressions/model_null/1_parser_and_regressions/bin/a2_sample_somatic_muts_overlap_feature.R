library(tidyverse)
library(data.table)
library(dtplyr)
library(GenomicRanges)
library(rtracklayer)
library(rlang)
library(valr)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("reduce", "IRanges")
conflict_prefer("expand", "tidyr")


### load sample; merge with dna repair, parse

# from command parameters
args = commandArgs(trailingOnly=TRUE)


# load utils.R (functions)
if(interactive()){
  source("/g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/bin/utils.R")
} else {
  source(args[1])
}

feature_name = ifelse(interactive(),
                 yes = "RepliSeq",
                 no = args[2])


sample_name = ifelse(interactive(),
                yes = "ABA07378_REP1", #AY9808_REP1"
                no = gsub("\\[|\\]", "", args[3])) # after channeling in nextflow, the sample names are contained within brackets, so remove them


somatic_data_paths = ifelse(interactive(),
                            yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/marcel_K562/2_rm_cups_lift_to_hg19/data_hg19/muts_pass_hg19_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/data/muts_pass_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/data/muts_pass_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/petljak_2022/parse_input/res/muts_pass_",
                            no = args[4]) %>% 
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1)


## load feature tracks file (already processed by first script)
feature_file_processed = ifelse(interactive(),
                                yes = Sys.glob("../work/*/*/feature_file_processed_RepliSeq.tsv")[1],
                                no = args[5]) %>%
  read_tsv %>% 
  select(seqnames, start, end, score) %>% 
  rename("chrom" = "seqnames")



## load sample (somatic mutations)

# find path+sample name that exists
existing_file = c()
for(file_exists in paste0(somatic_data_paths, sample_name, ".csv")){
  if(file.exists(file_exists)){
    existing_file = c(existing_file, file_exists)
  }
}

# raise error if no samples were found in any path, or if >=2 samples with same name exist in different paths (i.e. diff. datasets)
if(length(existing_file) != 1){
  stop(paste0("ERROR! No samples were found in any path OR multiple samples with the same name exist in different paths:", "\n", existing_file))
}


### load sample somatic muts, keep good mappability regions, merge to feature ranges, binarize score
sample_somatic_muts_overlap_feature = read_csv(existing_file) %>%
  select(chr, start, end, tri) %>% 
  rename("chrom" = "chr") %>%
  mutate(mut_id = paste0("mut_", row_number())) %>% 
  # merge to feature ranges, to get score
  bed_intersect(feature_file_processed, suffix = c("_SNVs", "_feature")) %>% 
  select(c(chrom, contains("_SNVs"), score_feature)) %>% 
  rename_all(~str_replace_all(., "_SNVs|_feature", "")) %>% 
  select(tri, score) %>% 
  table %>%
  as.data.frame
gc()

if(nrow(sample_somatic_muts_overlap_feature) == 0){
  sprintf("WARNING - Empty 'sample_somatic_muts_overlap_feature' table! Table will be empty...")
  sample_somatic_muts_overlap_feature = tibble("SBS96" = as.character(NA),
                                               "bin" = as.character(NA), 
                                               "mutcount" = 0,
                                               "sample_name" = sample_name)
} else {
  sample_somatic_muts_overlap_feature = sample_somatic_muts_overlap_feature %>%
    rename("SBS96" = "tri",
           "bin" = "score",
           "mutcount" = "Freq") %>% 
    arrange(SBS96, bin) %>% 
    mutate(sample_name = sample_name)
}

write_tsv(sample_somatic_muts_overlap_feature,
          paste0("sample_", sample_name, "_somatic_muts_overlap_", feature_name, ".tsv"))
