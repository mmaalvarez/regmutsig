library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("between", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("reduce", "IRanges")
conflict_prefer("desc", "dplyr")
conflict_prefer("reverseComplement", "spgs")
conflict_prefer("strsplit", "base")
conflict_prefer("expand", "tidyr")
conflict_prefer("first", "dplyr")
conflict_prefer("last", "dplyr")

# load utils.R (functions)
source("/g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/bin/utils.R")

# feature files paths
features_paths = read_csv("../../input_lists/features.csv")

dir.create("shuffled_feature_beds")

set.seed(1)

for(i in seq(1, nrow(features_paths))){
  
  feature_name_path = features_paths %>% 
    slice(i)
  feature_name = feature_name_path$name
  feature_path = feature_name_path$path

  # load feature tracks file
  feature_file = tryCatch(import.bw(feature_path),
                          error = function(e) tryCatch(import.bedGraph(feature_path),
                                                       error = function(e) tryCatch(import.bed(feature_path),
                                                                                    error = function(e) tryCatch(makeGRangesFromDataFrame(read_tsv(feature_path), keep.extra.columns = T),
                                                                                                                 error = function(e) import.wig(feature_path)))))
  
  names(elementMetadata(feature_file)) = "score"
  
  # shuffle score across ranges
  feature_file_shuffled = data.frame(feature_file) %>%
    transform(score = sample(score)) %>% 
    select(-c(width, strand))
  
  # collapse the newly with-same-score contiguous ranges
  feature_file_shuffled = unlist(reduce(split(makeGRangesFromDataFrame(feature_file_shuffled, keep.extra.columns = T), ~score)))
  # move back the scores to metadata column (they show as row names now) and rename as "score"
  mcols(feature_file_shuffled) = names(feature_file_shuffled)
  feature_file_shuffled = data.frame(feature_file_shuffled) %>%
    mutate(start = start-1) %>% 
    arrange(seqnames, start) %>%
    rename("score" = "X") %>% 
    select(-c(width, strand))
    
  write_tsv(feature_file_shuffled,
            paste0("shuffled_feature_beds/", feature_name, "_shuffled.bed"),
            col_names = F)
  gc()
}
