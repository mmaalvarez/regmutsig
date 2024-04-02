library(tidyverse)
library(msm) # rtnorm
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")


##### samples info

K562 = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/marcel_K562/metadata/WGS_clones_info.tsv") %>% 
  rename("Sample" = "sample_id",
         "alteration" = "genotype KOs") %>%
  mutate(dataset = "Marcel",
         Sample = gsub("_REP1", "", Sample),
         # I consider MSH3-/- as MMRwt
         altered_pathway_or_treatment_type = ifelse((str_detect(alteration, "MSH2|MSH6|MLH1|PMS2") | str_detect(`MMR deficiency expected`, "strong|middle")) & str_detect(alteration, "OGG1|MUTYH|NTH1|NEIL2"),
                                                    "MMR & BER",
                                                    ifelse(str_detect(alteration, "MSH2|MSH6|MLH1|PMS2") | str_detect(`MMR deficiency expected`, "strong|middle"),
                                                           "MMR",
                                                           ifelse(str_detect(alteration, "OGG1|MUTYH|NTH1|NEIL2"),
                                                                  "BER",
                                                                  "control"))),
         alteration = ifelse(alteration != "WT",
                             gsub("$", " KO", alteration),
                             alteration)) %>%  
  select(Sample, dataset, alteration, altered_pathway_or_treatment_type)

iPSC = c("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/sample_treatments.tsv",
         "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/sample_gene_ko.tsv") %>% 
  # only Sample and info* columns are selected
  map_df(~read_tsv(.x)) %>% 
  rename("alteration" = "info1",
         "altered_pathway_or_treatment_type" = "info2",
         "Sample" = "sample_id") %>% 
  mutate(altered_pathway_or_treatment_type = gsub("^[a-k]_", "", altered_pathway_or_treatment_type),
         altered_pathway_or_treatment_type = gsub("Control", "control", altered_pathway_or_treatment_type),
         dataset = ifelse(str_detect(Sample, "MSM0"),
                          "Kucab",
                          ifelse(str_detect(Sample, "MSK0"),
                                 "Zou",
                                 "ERROR: Unexpected sample name")),
         altered_pathway_or_treatment_type = gsub("DNA damage response inhibitors", "DNA damage resp. inh.", altered_pathway_or_treatment_type),
         # "MMRko" sample with very low mut burden, treat as control
         `altered_pathway_or_treatment_type` = ifelse(Sample == "MSK0.123_s1",
                                                      "control",
                                                      `altered_pathway_or_treatment_type`),
         alteration = ifelse(dataset == "Zou",
                             gsub("$", " KO", alteration),
                             alteration)) %>% 
  select(Sample, dataset, alteration, altered_pathway_or_treatment_type)

petljak = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/petljak_2022/info/metadata.tsv") %>% 
  rename("Sample" = "sample_id",
         "alteration" = "info1",
         "cell_line" = "info2") %>%
  mutate(dataset = "Petljak",
         # I consider MSH3-/- as MMRwt
         altered_pathway_or_treatment_type = ifelse(alteration == "WT",
                                                    "control",
                                                    "APOBEC"),
         alteration = gsub("_KO", " KO", alteration)) %>% 
  select(Sample, dataset, alteration, altered_pathway_or_treatment_type)

# merge datasets metadata
samples_info = bind_rows(K562, iPSC) %>% 
  bind_rows(petljak)



## load res

critical_value = 1.28 # for CI80%
# 1.96 for CI95%

coefficient_table = lapply(c("../../../1_parser_and_regressions/res/results_regional_feature_regressions_all_samples.tsv",
                             "../../../../model_2/1_parser_and_regressions/res/results_SBS96_regressions_all_samples.tsv"),
                           read_tsv) %>%
  Reduce(function(x, y) bind_rows(x, y), .) %>%
  # for regressions that died (estimate and/or CIs >=|10|), convert the estimate and its CI to 0, as these would mean that the regression died
  pivot_wider(names_from = stat, values_from = val) %>% 
  mutate(conf.low = estimate - (std.error*critical_value),
         conf.high = estimate + (std.error*critical_value),
         reg_died = ifelse(abs(estimate) >= 10 | conf.low <= -10 | conf.high >= 10, "yes", "no"),
         across(starts_with("estimate") | starts_with("conf."),
                ~ifelse(abs(estimate) >= 10 | conf.low <= -10 | conf.high >= 10, 0, .)),
         estimate = ifelse(is.na(estimate), 0, estimate),
         conf.low = ifelse(is.na(conf.low), 0, conf.low),
         conf.high = ifelse(is.na(conf.high), 0, conf.high)) %>%
  pivot_longer(cols = starts_with("estimate") | starts_with("conf."), names_to = "stat", values_to = "val") %>% 
  rename("Sample" = "sample_name") %>% 
  mutate(Sample = gsub("_REP1", "", Sample)) %>% 
  left_join(samples_info) %>% 
  select(Sample, dataset, alteration, altered_pathway_or_treatment_type, feature_name, non_ref_bin_feature_level, feature_type, stat, val) %>%
  pivot_wider(names_from = stat, values_from = val) %>% 
  group_by(Sample, feature_name)
gc()


## write original coefficients
coefficient_table %>% 
  select(Sample, feature_name, estimate) %>% 
  write_tsv("original_coeff.tsv")


## Generate matrices resampling betas from their CI % distributions

resample_from_CI = function(coefficient_table){
  summarise(coefficient_table, resampled_estimate = rtnorm(n = 1,
                                                           mean = estimate,
                                                           sd = 1, 
                                                           lower = conf.low,
                                                           upper = conf.high))
}

totalNumIters = 6000
coefficient_Resamp = list()
set.seed(1)

for (nIter in 1:totalNumIters) {
  print(paste0("nIter ", nIter))
  # for each sample (row) resample coefficients from CI % distrs.
  coefficient_Resamp[[nIter]] = resample_from_CI(coefficient_table) %>% 
    mutate("nIter" = nIter)
  gc()
}

# bind all permuted tables
coefficient_Resamp = bind_rows(coefficient_Resamp)

write_tsv(coefficient_Resamp, paste0("permuted_coefficients_", totalNumIters, "iters.tsv"))
