library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg19")
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")


args = commandArgs(trailingOnly=TRUE)


good_mappability_regions = ifelse(interactive(),
                                  yes = "/g/strcombio/fsupek_cancer3/malvarez/chromatin_info/good_mappability_regions/k50.umap.parsed.bed",
                                  no = args[1]) %>%
  import.bed() %>% data.frame %>%
  ## WARNING: here I keep 1 half of the genome
  filter(seqnames %in% paste0("chr", c(seq(1,21,2), "X"))) %>% 
  select(seqnames, start, end, width, strand) %>% 
  makeGRangesFromDataFrame()


# get sequence for each range
sequences = getSeq(BSgenome.Hsapiens.UCSC.hg19,
                   names = good_mappability_regions)
gc()


# get frequency of each trinucleotide per sequence, moving 1 nucleotide downstream each time
n_trinuc64_at_risk = trinucleotideFrequency(sequences) %>% 
  data.frame()
rm(sequences) ; gc()


#### "offset" is actually num trinuc32 At Risk (of each trinuc32) across all good mappability regions

AGT_column = c('>A)', '>G)', '>T)') %>% 
  data.frame %>% 
  `colnames<-`("AGT")

good_mappability_regions_n_trinuc32_at_risk = n_trinuc64_at_risk %>% 
  summarise_all(~sum(.)) %>% 
  ungroup %>% 
  pivot_longer(cols = everything(),
               names_to = 'trinuc64',
               values_to = "n_trinuc64_at_risk") %>%
  # trinucs that do not have a C or T in center, convert to reverse complement
  mutate(trinuc32 = ifelse(substr(trinuc64, start = 2, stop = 2) %in% c('A', 'G'),
                           spgs::reverseComplement(trinuc64, case="upper"),
                           trinuc64)) %>% 
  # sum up frequencies of each N(C|T)N & reverse complement pair, within id
  group_by(trinuc32) %>%
  summarise(n_trinuc32_at_risk = sum(n_trinuc64_at_risk)) %>%
  ungroup %>% 
  merge(AGT_column, all = T) %>% 
  mutate(SBS96 = paste0(substr(trinuc32, start = 1, stop = 1),
                        "(",
                        substr(trinuc32, start = 2, stop = 2),
                        AGT,
                        substr(trinuc32, start = 3, stop = 3)),
         # correct T>T to T>C
         SBS96 = gsub("T>T", "T>C", SBS96)) %>% 
  ungroup %>% 
  select(trinuc32, n_trinuc32_at_risk, SBS96) %>% 
  arrange(trinuc32, SBS96)
gc()

write_tsv(good_mappability_regions_n_trinuc32_at_risk,
          paste0("good_mappability_regions_n_trinuc32_at_risk.tsv"))
