library(tidyverse)
conflicts_prefer(dplyr::filter)

## Marcel K562 52 KO pairs (strong, mid, weak, or no dMMR, or BER)
K562_pairs = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/marcel_K562/metadata/WGS_clones_info.tsv") %>% 
  select(sample_id)

## Zou
zou_samples = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/sample_gene_ko.tsv") %>% 
  select(sample_id)

## Kucab
kucab_samples = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/sample_treatments.tsv") %>% 
  select(sample_id)

## Petljak
petljak_samples = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/petljak_2022/parse_input/input_list/sample_ids.tsv") %>% 
  select(sample_id)

## colorectal
colorectal = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/metadata/comb_metadata_final_6datasets__noconsent_samples_removed__hartwig_upd.tsv") %>%
	filter(OriginalType=="Colorectum" & tissue=="colorectum") %>%
	select(sample_id)


### bind all and write out
write_csv(bind_rows(bind_rows(bind_rows(bind_rows(K562_pairs, zou_samples), kucab_samples), petljak_samples), colorectal),
          "sample_names.csv", col_names = F)

