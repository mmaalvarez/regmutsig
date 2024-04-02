library(tidyverse)
library(ggrepel)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("Position", "ggplot2")
conflict_prefer("Position", "ggplot2")



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

coefficient_table = lapply(c("../../../../1_parser_and_regressions/res/results_regional_feature_regressions_all_samples.tsv",
                             "../../../../../model_2/1_parser_and_regressions/res/results_SBS96_regressions_all_samples.tsv"),
                           read_tsv) %>%
  Reduce(function(x, y) bind_rows(x, y), .) %>%
  # remove regional features that died (estimate and/or CIs >=|10|, it means that the regression died), to not be kept because they will have the lowest coefficient 
  pivot_wider(names_from = stat, values_from = val) %>% 
  mutate(conf.low = estimate - (std.error*critical_value),
         conf.high = estimate + (std.error*critical_value)) %>% 
  rowwise() %>% 
  mutate(reg_died = ifelse(abs(estimate) >= 10 | conf.low <= -10 | conf.high >= 10, "yes", "no")) %>%
  filter(reg_died == "no") %>% 
  pivot_longer(cols = starts_with("estimate") | starts_with("conf."), names_to = "stat", values_to = "val") %>% 
  rename("Sample" = "sample_name") %>% 
  mutate(Sample = gsub("_REP1", "", Sample)) %>% 
  left_join(samples_info) %>% 
  select(Sample, dataset, alteration, altered_pathway_or_treatment_type, feature_name, non_ref_bin_feature_level, feature_type, stat, val) %>%
  pivot_wider(names_from = stat, values_from = val)

write_tsv(coefficient_table, "coefficient_table.tsv")



## plot coefficient + CI % dists

jet_colors = colorRampPalette(c("red", "yellow", "green", "blue", "magenta", "black"))


## by dataset

for(dataset_name in unique(coefficient_table$dataset)){

  coefficient_table_boxplots = coefficient_table %>%
    # keep 1 dataset
    filter(dataset == dataset_name) %>%
    ## to plot only regional feature regressions
    filter(feature_type == "regional_feature") %>%
    # MBD4 has a lot of error
    filter(feature_name != "MBD4")

  axis_text_x_size = case_when(dataset_name %in% c("Zou", "Kucab") ~ 3, 
                               dataset_name %in% c("Marcel", "Petljak") ~ 6)
  
  plot_coefficient_CI = ggplot(coefficient_table_boxplots,
                               aes(x = feature_name,
                                   y = estimate)) +
    geom_errorbar(aes(ymin = `conf.low`,
                      ymax = `conf.high`),
                  linewidth = 0.1,
                  col = "red",
                  position = position_dodge2()) +
    geom_point(aes(col = feature_name),
               size = 0.3,
               alpha = 0.8,
               stat = "identity",
               position = position_dodge2(width = 0.9)) +
    scale_color_manual(values = jet_colors(length(unique(coefficient_table_boxplots$feature_name)))) +
    # show genotype of outlier samples
    geom_text_repel(data = coefficient_table_boxplots %>%
                      group_by(altered_pathway_or_treatment_type) %>%
                      slice_max(order_by = abs(conf.low) + conf.high,
                                n = 10,
                                with_ties = F),
              aes(label = alteration),
              size = 2,
              min.segment.length = 0.0000000000001,
              segment.size = 0.1,
              segment.linetype = "dashed",
              segment.color = "black",
              max.overlaps = 1000) +
    geom_hline(yintercept = 0, color = "darkgray", linetype = "dashed") +
    facet_wrap(facets = vars(altered_pathway_or_treatment_type),
               nrow = 4,
               scales = "free_y") +
    theme_classic() +
    xlab("Regional feature") +
    ylab("Regression coefficient + CI80%") +
    ggtitle("Samples altered pathway OR treatment type") +
    theme(text = element_text(size = 10),
          axis.text.x = element_text(size = axis_text_x_size, angle = 45, hjust = 1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 10))
  ggsave(paste0("plot_coefficient_CI_", dataset_name, ".jpg"),
         plot = plot_coefficient_CI,
         device = "jpg",
         width = 13,
         height = 7.3,
         dpi = 600)

  
  
  #############################################################
  ## plot individual coefficients in heatmap
  coefficient_table_heatmap = coefficient_table %>% 
    # keep 1 dataset
    filter(dataset == dataset_name) %>% 
    ## to plot only regional feature regressions
    filter(feature_type == "regional_feature") %>%
    # MBD4 has a lot of error
    filter(feature_name != "MBD4") %>% 
    # we only plot coeff, not CI
    select(-starts_with("conf.")) %>% 
    rename("term" = "feature_name",
           "Reg.\ncoeff." = "estimate") %>% 
    arrange(altered_pathway_or_treatment_type, alteration, term)
  
  coefficient_table_heatmap_plot = ggplot(coefficient_table_heatmap,
                                          aes(x = alteration,
                                              y = term)) +
    geom_tile(aes(fill = `Reg.\ncoeff.`)) +
    scale_fill_gradientn(colours = c('red', 'white', 'lightblue'),
                         values = scales::rescale(c(min(coefficient_table_heatmap$`Reg.\ncoeff.`), 
                                                    0, 
                                                    max(coefficient_table_heatmap$`Reg.\ncoeff.`)))) + 
    ## group sample labels in x axis by their altered_pathway_or_treatment_type
    facet_grid(cols = vars(altered_pathway_or_treatment_type),
               scales = "free", space = "free") +
    theme_bw() +
    xlab("Altered pathway OR treatment") +
    ylab("Regional feature") +
    theme(text = element_text(size = 10),
          axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 9),
          strip.text.x.top = element_text(size = 4),
          strip.background = element_blank())
  ggsave(paste0("coefficient_table_heatmap_", dataset_name, ".jpg"),
         plot = coefficient_table_heatmap_plot,
         device = "jpg",
         width = 13,
         height = 7,
         dpi = 600)
}
