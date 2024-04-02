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


## load res
coefficient_table = read_tsv("../../1_parser_and_regressions/res/results_regressions_full_simulations.tsv") %>%
  select(-theta) %>% 
  # remove regional features that died (estimate and/or CIs >=|10|, it means that the regression died), to not be kept because they will have the lowest coefficient
  rowwise() %>% 
  filter(!(abs(estimate) >= 10 | conf.low <= -10 | conf.high >= 10))

phenotype_levels = c("wild-type", paste0("def-feature",seq(1,length(unique(coefficient_table$sample_type))-1)))
feature_levels = paste0("feature", seq(1,length(unique(coefficient_table$term))))

coefficient_table = coefficient_table %>% 
  mutate(phenotype = factor(sample_type, levels = phenotype_levels, ordered = T),
         term = factor(term, levels = feature_levels, ordered = T)) %>% 
  select(parameters, sample_id, phenotype, term, estimate, starts_with("conf.")) %>% 
  separate(parameters, into = c('nFeatures',
                                'nBins',
                                'binsize',
                                'nDefRepSamplesPerFeature',
                                'nWTsamples',
                                'mutrate_fully_unrepaired_bins',
                                'mutrate_fully_repaired_bins'), sep = "_") %>% 
  # mutrates as scientific notation
  mutate(mutrate_fully_unrepaired_bins = formatC(as.numeric(mutrate_fully_unrepaired_bins), format = "e", digits = 0),
         mutrate_fully_repaired_bins = formatC(as.numeric(mutrate_fully_repaired_bins), format = "e", digits = 0))
write_tsv(coefficient_table, "coefficient_table.tsv")


#############################################################
## plot individual coefficients in heatmap

# select subset of some parameters
coefficient_table_subset = coefficient_table %>% 
  filter(nFeatures %in% c(2,10)) %>% 
  filter(nDefRepSamplesPerFeature %in% c(1,20)) %>%
  filter(nWTsamples %in% c(25,200))

coefficient_table_heatmap_plot = ggplot(coefficient_table_subset,
                                        aes(x = phenotype,
                                            y = term)) +
  geom_tile(aes(fill = estimate)) +
  scale_fill_gradientn(colours = c('red', 'white', 'lightblue'),
                       values = scales::rescale(c(min(coefficient_table_subset$estimate), 
                                                  0, 
                                                  max(coefficient_table_subset$estimate))),
                       name = "Regression coefficient") + 
  ggh4x::facet_nested(cols = vars(mutrate_fully_unrepaired_bins,
                                  mutrate_fully_repaired_bins),
                      rows = vars(nFeatures,
                                  nDefRepSamplesPerFeature,
                                  nWTsamples),
                      space = "free", scales = "free", #independent = "x",
                      labeller = label_both) +
  theme_bw() +
  xlab("Sample phenotype") +
  ylab("Regression variable ('low' vs. 'high' regional activity)") +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 5),
        legend.position = "top")
ggsave("coefficient_table_heatmap.jpg",
       plot = coefficient_table_heatmap_plot,
       device = "jpg",
       width = 23,
       height = 13,
       dpi = 600)


####################################################
## plot coefficient + CI % dists

## IT IS TOO MUCH for it to compute, maybe should remove CIs, and/or dodging, and/or facets

# jet_colors = colorRampPalette(c("red", "yellow", "green", "blue", "magenta", "black"))

# plot_coefficient_CI = ggplot(coefficient_table,
#                              aes(x = term,
#                                  y = estimate)) +
#   geom_hline(yintercept = 0, color = "darkgray", linetype = "dashed") +
#   geom_errorbar(aes(ymin = `conf.low`, 
#                     ymax = `conf.high`), 
#                 linewidth = 0.1,
#                 col = "darkgray",
#                 position = position_dodge2()) +
#   geom_point(aes(col = term),
#              size = 0.1,
#              stat = "identity",
#              position = position_dodge2(width = 0.9)) +
#   scale_color_manual(values = jet_colors(length(unique(coefficient_table$term)))) +
#   # show phenotype of outlier samples
#   geom_text_repel(data = coefficient_table %>% 
#                     filter(conf.low < -2 | conf.high > 2) %>%
#                     group_by(term),
#                   aes(label = phenotype),
#                   size = 1,
#                   max.overlaps = 1000,
#                   min.segment.length = 10000000) +
#   facet_wrap(facets = vars(phenotype), nrow = 1) +
#   # ggh4x::facet_nested(rows = vars(mutrate_fully_unrepaired_bins,
#   #                                  mutrate_fully_repaired_bins,
#   #                                  phenotype),
#   #                        cols = vars(number_features,
#   #                                    number_deficient_repair_samples_per_feature,
#   #                                    total_number_wildtype_samples),
#   #                        space = "free",
#   #                        scales = "free") +
#   theme_bw() +
#   xlab("Regression variable") +
#   ylab("Regression coefficient + CI80%") +
#   theme(text = element_text(size = 10),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.grid.minor.y = element_blank(),
#         legend.position = "none")
# ggsave("plot_coefficient_CI.jpg",
#        plot = plot_coefficient_CI,
#        device = "jpg",
#        width = 13,
#        height = 7.3,
#        dpi = 600)