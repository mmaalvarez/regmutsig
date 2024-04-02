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
                      yes = "MBD4",
                      no = args[2])

## load feature tracks file with good mappability
feature_file_good_mappability = ifelse(interactive(),
                                       yes = Sys.glob(paste0("../work/*/*/",feature_name,"_good_mappability.tsv"))[1],
                                       no = args[3]) %>% 
  read_tsv

## load score used for cutoff of high (vs low)
cutoff_score = ifelse(interactive(),
                      yes = Sys.glob(paste0("../work/*/*/cutoff_score_",feature_name,".tsv"))[1],
                      no = args[4]) %>% 
  read_tsv %>% 
  pull(cutoff_score)

## to keep SNVs in good mappability regions
good_mappability_regions = ifelse(interactive(),
                                  yes = "/g/strcombio/fsupek_cancer3/malvarez/chromatin_info/good_mappability_regions/k50.umap.parsed.bed",
                                  no = args[5]) %>%
  import.bed() %>% data.frame %>%
  select(seqnames, start, end, width, strand) %>% 
  makeGRangesFromDataFrame()

## parameter values for trinuc_matching()
stoppingCriterion = ifelse(interactive(),
                           yes = "0.01",
                           no = args[6]) %>% 
  as.numeric()

max_fraction_removed_trinucs = ifelse(interactive(),
                                      yes = "0.2",
                                      no = args[7]) %>% 
  as.numeric()

# bp 
flank_length = ifelse(interactive(),
                      yes = "50000",
                      no = args[8]) %>% 
  as.numeric()



# apply cutoff for the feature of this iteration, for binarization (and later collapse into 2 bins) by being lower or larger than the cutoff
if(is.numeric(cutoff_score)){
  
  feature_file_good_mappability_binarized = feature_file_good_mappability %>% 
    lazy_dt %>% 
    mutate(score = ifelse(score >= cutoff_score,
                          "high", "low")) %>% 
    as_tibble %>% 
    makeGRangesFromDataFrame(keep.extra.columns = T)
  
} else { # scores are already categories
  
  feature_file_good_mappability_binarized = feature_file_good_mappability %>% 
    # both score levels to 'low' and 'high', in case they are otherwise
    mutate(score = str_replace_all(score, 
                                   c("0\\-3" = "low", 
                                     "0\\&1" = "low", 
                                     "0" = "low", 
                                     "bgGenome" = "low",
                                     "4\\-6" = "high", 
                                     "2\\&3" = "high", 
                                     "1\\-3" = "high"))) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = T)
}
rm(feature_file_good_mappability) ; gc()


## collapse contiguous ranges with same level
feature_file_good_mappability_binarized_collapsed = unlist(reduce(split(feature_file_good_mappability_binarized, ~score)))
# move back the scores to metadata column (they show as row names now) and rename as "score"
mcols(feature_file_good_mappability_binarized_collapsed) = names(feature_file_good_mappability_binarized_collapsed)
feature_file_good_mappability_binarized_collapsed = data.frame(feature_file_good_mappability_binarized_collapsed) %>% 
  arrange(seqnames, start) %>%
  rename("score" = "X")
rm(feature_file_good_mappability_binarized) ; gc()



### NEW: restrict the "low" feature level genome to just the flank_length bp flanking each "high" region
# if a low row is shorter than 2*flank_length, kept as is
# afterwards the low-mappability regions within these flanks will be removed again (as they were readded to simplify the process)

feature_file_good_mappability_binarized_collapsed_withflanks = feature_file_good_mappability_binarized_collapsed %>%
  # fully collapse contiguous "low" rows, and keep other rows as they are
  mutate(contiguous_low = cumsum(score != lag(score, default = first(score))) * (score == "low") + 1,
         # "high" rows will have a 0, so they are ignored in the next step
         contiguous_low = ifelse(score == "high",
                                 0,
                                 contiguous_low)) %>% 
  group_by(seqnames, contiguous_low) %>%
  mutate(start = ifelse(contiguous_low != 0,
                        first(start),
                        start),
         end = ifelse(contiguous_low != 0,
                      last(end),
                      end),
         score = ifelse(contiguous_low != 0,
                        first(score),
                        score)) %>%
  ungroup %>% 
  select(-c(width, contiguous_low)) %>% 
  distinct %>% 
  mutate(width = end-start+1) %>% 
  relocate(width, .before = strand) %>% 
  arrange(seqnames, start) %>% 
  # label rows by type
  mutate(type_row = ifelse(row_number() == 1,
                           "first row",
                           ifelse(row_number() == nrow(.),
                                  "last row",
                                  ifelse(score == "high",
                                         "high score",
                                         ifelse(width < flank_length*2,
                                                "short flanks",
                                                # these are the ones that we will trim:
                                                "long flanks")))),
         # deal with possible low rows in the first and last positions
         start = ifelse(type_row == "first row" & score == "low" & width > flank_length,
                        end - flank_length + 1,
                        start),
         end = ifelse(type_row == "last row" & score == "low" & width > flank_length,
                      start + flank_length - 1,
                      end),
         width = end-start+1)
rm(feature_file_good_mappability_binarized_collapsed) ; gc()

# duplicate and label "long flanks" rows
tmp = feature_file_good_mappability_binarized_collapsed_withflanks %>%
  filter(type_row == "long flanks") %>%
  # Use 'uncount' to duplicate rows
  uncount(2, .id = "dup_id") %>%
  unite("type_row", type_row, dup_id, sep = " ")

# readd to the rest of the data
feature_file_good_mappability_binarized_collapsed_withflanks = feature_file_good_mappability_binarized_collapsed_withflanks %>%
  select(-type_row) %>% 
  left_join(tmp) %>% 
  # trim flanks
  mutate(start = ifelse(!is.na(type_row) & type_row == "long flanks 2",
                        end - flank_length + 1,
                        start),
         end = ifelse(!is.na(type_row) & type_row == "long flanks 1",
                      start + flank_length - 1,
                      end),
         width = end-start+1) %>% 
  select(-type_row)
rm(tmp) ; gc()


## again, keep only good_mappability regions from the "low" rows
feature_file_good_mappability_binarized_collapsed_withflanks_list = list()

existing_normal_chromosomes = unique(feature_file_good_mappability_binarized_collapsed_withflanks$seqnames)[unique(feature_file_good_mappability_binarized_collapsed_withflanks$seqnames) %in% paste0("chr", c(seq(1,22), "X", "Y"))]

for(chr in existing_normal_chromosomes){
  
  # split high and low rows per chromosome, as the high rows are already good mappability only
  single_chr_high_rows = feature_file_good_mappability_binarized_collapsed_withflanks %>% 
    filter(seqnames == chr) %>% 
    filter(score == "high") %>% 
    ## add 1 bp downstream and upstream regions of width 1 or 2, so trinucs can be fetched later
    mutate(width = end-start+1,
           start = ifelse(width <= 2,
                          start-1,
                          start),
           # if width was 1, after adding 1bp to start it will still be 2bp, so add another bp to end
           end = ifelse(width == 1,
                        end+1,
                        end),
           width = end-start+1,
           strand = "*") %>% 
    relocate(score, .after = strand)
  
  single_chr_low_rows = feature_file_good_mappability_binarized_collapsed_withflanks %>% 
    filter(seqnames == chr) %>% 
    filter(score == "low") %>% 
    makeGRangesFromDataFrame(keep.extra.columns = T)
  
  single_chr_low_rows_good_mappability = findOverlaps_wrapper(single_chr_low_rows,
                                                              good_mappability_regions) %>%
    ## add 1 bp downstream and upstream regions of width 1 or 2, so trinucs can be fetched later
    mutate(width = end-start+1,
           start = ifelse(width <= 2,
                          start-1,
                          start),
           # if width was 1, after adding 1bp to start it will still be 2bp, so add another bp to end
           end = ifelse(width == 1,
                        end+1,
                        end),
           width = end-start+1,
           strand = "*") %>% 
    relocate(score, .after = strand)
  
  feature_file_good_mappability_binarized_collapsed_withflanks_list[[chr]] = bind_rows(single_chr_high_rows, single_chr_low_rows_good_mappability) %>% 
    arrange(seqnames, start)
  
  rm(single_chr_high_rows) ; rm(single_chr_low_rows) ; rm(single_chr_low_rows_good_mappability) ; gc()
}

feature_file_final = bind_rows(feature_file_good_mappability_binarized_collapsed_withflanks_list)
rm(feature_file_good_mappability_binarized_collapsed_withflanks) ; rm(feature_file_good_mappability_binarized_collapsed_withflanks_list) ; rm(good_mappability_regions) ; gc()

write_tsv(feature_file_final, paste0("feature_file_processed_", feature_name,".tsv"))



## get sequence for each range
sequences = getSeq(BSgenome.Hsapiens.UCSC.hg19,
                   names = makeGRangesFromDataFrame(feature_file_final, 
                                                    keep.extra.columns = T))
gc()

# get frequency of each trinucleotide per sequence, moving 1 nucleotide downstream each time
trinuc64_freq = trinucleotideFrequency(sequences) %>% 
  data.frame()
rm(sequences) ; gc()


trinuc_matching_input = feature_file_final %>% 
  as_tibble %>% 
  select(score) %>% 
  ## bind trinuc64_freq
  bind_cols(trinuc64_freq) %>% 
  group_by(score) %>% 
  summarise_all(~sum(.)) %>% 
  ungroup %>% 
  pivot_longer(cols = -c(score),
               names_to = 'trinuc64',
               values_to = "freq") %>%
  group_by(score) %>%
  # trinucs that do not have a C or T in center, convert to reverse complement
  mutate(trinuc32 = ifelse(substr(trinuc64, start = 2, stop = 2) %in% c('A', 'G'),
                           spgs::reverseComplement(trinuc64, case="upper"),
                           trinuc64)) %>% 
  # sum up frequencies of each N(C|T)N & reverse complement pair, within id
  group_by(score, trinuc32) %>%
  summarise(freq = sum(freq)) %>%
  ungroup %>% 
  # back to original format (each row ('id') maps to the same row in map_features)
  pivot_wider(names_from = 'trinuc32',
              values_from = 'freq') %>% 
  arrange(score) %>% 
  column_to_rownames("score")
rm(feature_file_final) ; rm(trinuc64_freq) ; gc()

bin_names = data.frame("bin" = rownames(trinuc_matching_input)) 


## run matching
matched_trinuc_counts = trinuc_matching(trinuc_matching_input, 
                                        # desired Euclidean score (max. overall distance between any bin's trinuc frequencies and all-bin-average trinuc frequencies; default 0.01)
                                        stoppingCriterion = stoppingCriterion,
                                        # don't allow to remove more total trinucleotide counts than this fraction of the total original trinucleotide counts (default 0.2)
                                        max_fraction_removed_trinucs = max_fraction_removed_trinucs) 
gc()

euclidean_score = matched_trinuc_counts[2][[1]]

matched_trinuc_counts = matched_trinuc_counts[1][[1]]


#### "offset" is actually num trinuc At Risk (per bin)

## sum the remaining counts per bin (high/low), so it's input for offset
offset_temp = matched_trinuc_counts %>%
  # sum up trinuc32 frequencies within bin profile
  column_to_rownames("bin") %>% 
  mutate(trinucAtRisk = rowSums(.)) %>% 
  rownames_to_column("bin") %>% 
  pivot_longer(cols = matches("^[A,C,G,T][C,T][A,C,G,T]$"),
               names_to = 'trinuc32',
               values_to = "n_trinuc32_at_risk_in_bin") %>% 
  # triplicate each row, adding '( >A)', '( >G)' and '( >T)' around the central C or T
  group_by_at(vars(!matches("n_trinuc32_at_risk_in_bin"))) %>% 
  slice(rep(row_number(), 3))

AGT_column = rep(c('>A)', '>G)', '>T)'), 
                 times = length(rownames(offset_temp)) / 3) %>% 
  data.frame %>% 
  `colnames<-`("AGT")

offset_temp = offset_temp %>% 
  bind_cols(AGT_column) %>% 
  mutate(SBS96 = paste0(substr(trinuc32, start = 1, stop = 1),
                        "(",
                        substr(trinuc32, start = 2, stop = 2),
                        AGT,
                        substr(trinuc32, start = 3, stop = 3)),
         # correct T>T to T>C
         SBS96 = gsub("T>T", "T>C", SBS96)) %>% 
  ungroup %>% 
  select(-AGT)
gc()

### we also want the fraction of each trinuc type that was removed per bin (high/low), for removing randomly from somatic mutations
prematching_trinuc32_freqs = trinuc_matching_input %>% 
  rownames_to_column("bin") %>% 
  pivot_longer(cols = matches("^[A,C,G,T][C,T][A,C,G,T]$"),
               names_to = 'trinuc32',
               values_to = "PREMATCHING_n_trinuc32_removed_in_bin")

offset = trinuc_matching_input %>%
  # subtract orig from remaining counts
  map2_dfc(column_to_rownames(matched_trinuc_counts, "bin"), ~ .x - .y) %>% 
  # bin names back, and trinucAtRisk
  bind_cols(bin_names) %>% 
  relocate(bin) %>% 
  pivot_longer(cols = matches("^[A,C,G,T][C,T][A,C,G,T]$"),
               names_to = 'trinuc32',
               values_to = "n_trinuc32_removed_in_bin") %>% 
  left_join(offset_temp) %>% 
  left_join(prematching_trinuc32_freqs) %>% 
  mutate(fraction_trinuc32_removed_in_bin = n_trinuc32_removed_in_bin / PREMATCHING_n_trinuc32_removed_in_bin,
         feature_name = feature_name,
         euclidean_score = euclidean_score) %>% 
  select(bin, trinucAtRisk, trinuc32, n_trinuc32_at_risk_in_bin, PREMATCHING_n_trinuc32_removed_in_bin, fraction_trinuc32_removed_in_bin, SBS96, feature_name, euclidean_score)
gc()

write_tsv(offset,
          paste0("trinuc_fractions_rm_per_bin_and_trinucAtRisk_", feature_name, ".tsv"))
