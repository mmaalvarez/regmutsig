library(tidyverse)
library(data.table)
library(dtplyr)
library(GenomicRanges)
library(rtracklayer)
library(valr)
library("BSgenome.Hsapiens.UCSC.hg19")
library(spgs)
library(rlang)
library(Biostrings)
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


args = commandArgs(trailingOnly=TRUE)

# load utils.R (functions)
if(interactive()){
  source("/g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/bin/utils.R")
} else {
  source(args[1])
}

feature_name = ifelse(interactive(),
                      yes = "SETD2",
                      no = args[2])

feature_path = ifelse(interactive(),
                      yes = "/g/strcombio/fsupek_cancer3/malvarez/chromatin_info/protein_binding/DNA_repair/MMR/MSH6_SETD2_H3K36me3_guominli/2_fold_enrichment_vs_input/bedgraphs_fold_enrich/collapse_ranges/SETD2/SRR13314134_sorted_FE.all_chr.collapsed_ranges.bed.gz",
                      no = args[3]) %>% 
  # remove comment
  gsub("( |\t).*", "", .)

## to keep SNVs in good mappability regions
good_mappability_regions = ifelse(interactive(),
                                  yes = "/g/strcombio/fsupek_cancer3/malvarez/chromatin_info/good_mappability_regions/k50.umap.parsed.bed",
                                  no = args[4]) %>%
  import.bed() %>% data.frame %>%
  select(seqnames, start, end, width, strand) %>% 
  makeGRangesFromDataFrame()


## list of ALL samples (to calculate feature's low-high score cutoff)
sample_names_list = ifelse(interactive(),
                           yes = "../input_lists/sample_names.csv",
                           no = args[5]) %>% 
  read_csv(col_names = F) %>% pull(X1)

somatic_data_paths = ifelse(interactive(),
                            yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/marcel_K562/2_rm_cups_lift_to_hg19/data_hg19/muts_pass_hg19_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/data/muts_pass_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/data/muts_pass_,/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/petljak_2022/parse_input/res/muts_pass_",
                            no = args[6]) %>% 
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1)

## earliest RT (level 6) in order to remove the potential effect of RT on the decrease of mut rate as a feature score increases (since many features' "high" scores could be correlated with early RT)
earliest_RT = ifelse(interactive(),
                     yes = "/g/strcombio/fsupek_data/_LEGACY/epi/replicationTiming/percentileBeds2016/RepliSeq.mean8solid__eqFreqBin6of6_avgSm0_lowThr0.00.bed.gz",
                     no = args[7]) %>%
  import.bed() %>% data.frame %>%
  select(seqnames, start, end, width, strand) %>% 
  makeGRangesFromDataFrame()


## load feature tracks file
feature_file = tryCatch(import.bw(feature_path),
                        error = function(e) tryCatch(import.bedGraph(feature_path),
                                                     error = function(e) tryCatch(import.bed(feature_path),
                                                                                  error = function(e) tryCatch(makeGRangesFromDataFrame(read_tsv(feature_path), keep.extra.columns = T),
                                                                                                               error = function(e) import.wig(feature_path)))))

names(elementMetadata(feature_file)) = "score"

feature_file = data.frame(feature_file) #%>% filter(seqnames %in% c("chr21", "chr22"))
scores = feature_file$score

## if a) "score" is numeric... OR b) if it only contains numbers, even though it is labelled as a "character" column and these numbers are formatted as strings
if(is.numeric(scores)       |       all(grepl("^[-]?[0-9]+(\\.[0-9]+)?$", unique(scores)))){
  
  ## ...and if max(scores)<1, scale between 0 and 1, e.g. for RMDflat...
  if(max(scores)<1){
    scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
    feature_file$score = scale_values(feature_file$score)
  }

  ## ...and/or if scores are too continuous, round to 1 decimal and collapse 
  if(length(unique(feature_file$score)) > 1000){
    feature_file$score = round(feature_file$score, 1)
    feature_file = unlist(reduce(split(makeGRangesFromDataFrame(feature_file, keep.extra.columns = T), ~score)))
    gc()
    # move back the scores to metadata column (they show as row names now) and rename as "score"
    mcols(feature_file) = names(feature_file)
    feature_file = data.frame(feature_file) %>%
      arrange(seqnames, start) %>%
      rename("score" = "X") %>%
      mutate(score = as.numeric(score))
  }
}


### keep good mappability regions in feature

# loop through chromosomes, otherwise it requires too much memory
feature_file_good_mappability_list = list()
feature_file_good_mappability_earliestRT_list = list()

existing_normal_chromosomes = unique(feature_file$seqnames)[unique(feature_file$seqnames) %in% paste0("chr", c(seq(1,22), "X", "Y"))]

for(chr in existing_normal_chromosomes){
  
  feature_file_single_chr = feature_file %>% 
    filter(seqnames == chr) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = T)
  
  feature_file_single_chr_good_mappability = findOverlaps_wrapper(feature_file_single_chr, # %>% data.frame %>% filter(seqnames=="chr21") %>% makeGRangesFromDataFrame(keep.extra.columns = T),   
                                                                  good_mappability_regions) %>% 
    mutate(width = end-start+1)
  feature_file_good_mappability_list[[chr]] = feature_file_single_chr_good_mappability
  
  # only extract the repliseq level 6 IF it's a numeric score feature, otherwise it's not necessary because we don't need to define a cutoff score value
  if(is.numeric(scores)       |       all(grepl("^[-]?[0-9]+(\\.[0-9]+)?$", unique(scores)))){
    feature_file_single_chr_good_mappability_earliestRT = findOverlaps_wrapper(makeGRangesFromDataFrame(feature_file_single_chr_good_mappability, 
                                                                                                        keep.extra.columns = T), # %>% data.frame %>% filter(seqnames=="chr21") %>% makeGRangesFromDataFrame(keep.extra.columns = T),   
                                                                               earliest_RT) %>% 
      mutate(width = end-start+1)
    feature_file_good_mappability_earliestRT_list[[chr]] = feature_file_single_chr_good_mappability_earliestRT
    rm(feature_file_single_chr_good_mappability_earliestRT)
  }
  rm(feature_file_single_chr) ; rm(feature_file_single_chr_good_mappability) ; gc()
}

# for next processes
feature_file_good_mappability = bind_rows(feature_file_good_mappability_list)
write_tsv(feature_file_good_mappability,
          paste0(feature_name, "_good_mappability.tsv"))

rm(feature_file) ; rm(feature_file_good_mappability_list) ; gc()


# AGAIN, if a) "score" is numeric... OR b) if it only contains numbers, even though it is labelled as a "character" column and these numbers are formatted as strings

scores = feature_file_good_mappability$score

if(is.numeric(scores)       |       all(grepl("^[-]?[0-9]+(\\.[0-9]+)?$", unique(scores)))){

  feature_file_good_mappability_earliestRT = bind_rows(feature_file_good_mappability_earliestRT_list)
  rm(feature_file_good_mappability_earliestRT_list)
  
  ## trim tail of scores (>quantile_val% quantile)
  quantile_val = 0.95
  score_q = quantile(scores, probs = quantile_val)
  feature_file_good_mappability_earliestRT = feature_file_good_mappability_earliestRT %>% 
    filter(score <= score_q)
  
  ## find path + sample names that exist
  samples_mutburdens = tibble("score" = double(),
                              "mutcount" = integer(),
                              "Sample" = character())
  
  for(sample_name in sample_names_list){
    
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
    
    sample_somatic_muts = read_csv(existing_file) %>%
      select(chr, start, end, tri) %>% 
      rename("chrom" = "chr") %>%
      mutate(mut_id = paste0("mut_", row_number())) %>% 
      # merge to feature ranges, to map SNVs to scores
      bed_intersect(rename(feature_file_good_mappability_earliestRT,
                           "chrom" = "seqnames"),
                    suffix = c("_SNVs", "_feature")) %>% 
      select(c(chrom, contains("_SNVs"), score_feature)) %>% 
      rename_all(~str_replace_all(., "_SNVs|_feature", "")) %>% 
      select(score) %>% 
      table %>%
      as.data.frame

    ## in case there are no mutations left after mapping to earliestRT
    if(nrow(sample_somatic_muts)!=0){
      
      sample_somatic_muts = sample_somatic_muts %>% 
        `colnames<-`(c("score", "mutcount")) %>% 
        mutate(Sample = sample_name,
               score = as.double(as.character((score))))
    
      samples_mutburdens = bind_rows(samples_mutburdens,
                                     sample_somatic_muts)
    }
    gc()
  }
  
  ## normalize mutcounts dividing by a given sample's total mut. burden and then by the total bp of a given score bin
  samples_mutburdens_total = samples_mutburdens %>% 
    group_by(Sample) %>% 
    summarise(sample_mutburden = sum(mutcount))
  
  score_bin_lengths_bp_log = feature_file_good_mappability_earliestRT %>% 
    group_by(score) %>% 
    summarise(log_score_bin_bp = log(sum(width)) + 1) %>% 
    mutate(score = as.character(round(score,4)))
  
  samples_mutburdens_normalized = merge(samples_mutburdens,
                                        samples_mutburdens_total) %>% 
    mutate(score = as.character(round(score,4))) %>% 
    merge(score_bin_lengths_bp_log) %>% 
    rowwise() %>% 
    mutate(score = as.numeric(score),
           mutcount_normalized = mutcount / sample_mutburden / log_score_bin_bp)
  
  write_tsv(samples_mutburdens_normalized,
            paste0(feature_name, "_scorewise_samples_mutburdens_normalized.tsv"))
  #samples_mutburdens_normalized = read_tsv(Sys.glob(paste0("../work/*/*/", feature_name, "_scorewise_samples_mutburdens_normalized.tsv"))[1])
  
  # function to find the "elbow", uses the central difference to approximate the 2nd derivative -- stackoverflow.com/a/4473065
  centralDiff = function(data){
    
    mean_mutcounts = data %>% 
      group_by(score) %>% 
      summarise(mean_mutcount = mean(mutcount_normalized)) %>% 
      arrange(score)
    
    scores_vector = mean_mutcounts %>% 
      pull(score)
    mean_mutcounts_vector = mean_mutcounts %>% 
      pull(mean_mutcount)
    
    central_differences = c()
    
    for(i in seq(which.max(mean_mutcounts_vector), length(mean_mutcounts_vector))){
      if(i==which.max(mean_mutcounts_vector)){
        central_difference_i = mean_mutcounts_vector[i+1] - mean_mutcounts_vector[i]
      } else if(i==length(mean_mutcounts_vector)){
        central_difference_i = 1
      } else{
        central_difference_i = mean_mutcounts_vector[i+1] - mean_mutcounts_vector[i-1] - mean_mutcounts_vector[i]
      }
      central_differences = c(central_differences, central_difference_i)
    }
    
    i_min = which.min(central_differences)
    
    # genomic ranges with a score equal or larger than score_high_cutoff will be classified as "high"
    score_high_cutoff = scores_vector[i_min + which.max(mean_mutcounts_vector)]
  }
  
  score_high_cutoff = centralDiff(samples_mutburdens_normalized)
  
  ## get low vs. high score cutoff ("elbow" -2nd derivative- in mutcounts~score)
  mutburdens_normalized_plot = ggplot(samples_mutburdens_normalized,
                                      aes(x = score,
                                          y = mutcount_normalized)) +
    scale_x_continuous(breaks = seq(0, round(max(scores), 0), 0.5)) +
    geom_boxplot(aes(group = score)) +
    geom_violin(aes(group = score)) +
    geom_point() +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cr")) +
    geom_vline(xintercept = score_high_cutoff, linetype = "dashed", color = "red") +
    theme_bw() +
    xlab(paste0(feature_name, " FE score (cap at ", quantile_val*100, "th quantile)")) +
    ylab("Sample's #SNVs / sample's mut. burden / log(score bin's length)")
  ggsave(paste0(feature_name, ".score_vs_norm.mut.burdens.jpg"),
         plot = mutburdens_normalized_plot,
         device = "jpg",
         width = 10,
         height = 5.6,
         dpi = 600)
  
} else { # if non-numerical levels, kept as they are
  score_high_cutoff = paste(sort(unique(scores)), collapse = ",")
  # create dummy empty plot
  write(data.frame(),
        paste0("EMPTY_categorical_", feature_name, ".score_vs_norm.mut.burdens.jpg"))
}

write_tsv(data.frame(feature_name, 
                     "cutoff_score" = score_high_cutoff),
          paste0("cutoff_score_",feature_name,".tsv"))
