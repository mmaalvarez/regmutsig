library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("Position", "ggplot2")

# ### THEORETICAL DISTRIBUTIONS
# n = 100000
# 
# ideal_normal_sample = rnorm(n, mean=0, sd=1)
# png("1_ideal_normal_sample.png")
# hist(ideal_normal_sample, breaks = 250)
# dev.off()
# 
# ideal_uniform_sample = runif(n, min=0, max=1)
# png("2_ideal_uniform_sample.png")
# hist(ideal_uniform_sample, breaks = 250)
# dev.off()
# 
# merged_ideal_samples = data.frame("pval" = ideal_normal_sample) %>% 
#   mutate(type = "ideal normal",
#          feature_name = "ideal normal") %>% 
#   bind_rows(data.frame("pval" = ideal_uniform_sample) %>% 
#               mutate(type = "ideal uniform",
#                      feature_name = "ideal uniform"))


#####################################################
### ACTUAL ANALYSIS

### load null (shuffled-feature) pvalues
null_pvals = read_tsv("../../../../1_parser_and_regressions/model_null/1_parser_and_regressions/res/results_regional_feature_regressions_all_samples.tsv") %>% 
  filter(stat == "p.value") %>% 
  rename("pval" = "val") %>% 
  mutate(type = "shuffled features") %>% 
  select(pval, type, feature_name)
png("3_null_pvals.png")
hist(null_pvals$pval, breaks = 250)
dev.off()

### load real pvalues
real_pvals = read_tsv("../../../../1_parser_and_regressions/res/results_regional_feature_regressions_all_samples.tsv") %>% 
  filter(stat == "p.value") %>% 
  rename("pval" = "val") %>% 
  mutate(type = "real features") %>% 
  select(pval, type, feature_name)
png("4_real_pvals.png")
hist(real_pvals$pval, breaks = 250)
dev.off()


# ### qqplots 
# 
# all_pvals = bind_rows(null_pvals,
#                       real_pvals,
#                       merged_ideal_samples)
#
# ## the sample from the ideal uniform AND the shuffled features pvalues should fit the qqplot
# ## this is the case also for the real features that DON'T have signal
# qqplot_uniform = ggplot(all_pvals,
#                         aes(sample = pval,
#                             color = type)) +
#   geom_qq(distribution = qunif,
#                          alpha = 0.5,
#                          size = 0.5,
#                          shape = 4) +
#   scale_color_manual(values = c("darkgreen", "yellow", "red", "blue")) +
#   guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   facet_wrap(facets = vars(feature_name)) +
#   xlab("Theoretical (uniform)") +
#   theme_classic() +
#   theme(legend.title = element_blank(),
#         legend.position = "top")
# ggsave("5_qqplot_uniform.jpg",
#        plot = qqplot_uniform,
#        device = "jpg",
#        width = 10,
#        height = 5.6,
#        dpi = 600,
#        bg = "white")
# 
# ## in the case of the real features that actually have some signal, their pvalues should behave more like the ideal normal
# qqplot_normal = ggplot(all_pvals,
#                        aes(sample = pval,
#                            color = type)) +
#   geom_qq(distribution = qnorm,
#           alpha = 0.5,
#           size = 0.5,
#           shape = 4) +
#   scale_color_manual(values = c("darkgreen", "yellow", "red", "blue")) +
#   guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   facet_wrap(facets = vars(feature_name)) +
#   xlab("Theoretical (normal)") +
#   theme_classic() +
#   theme(legend.title = element_blank(),
#         legend.position = "top")
# ggsave("6_qqplot_normal.jpg",
#        plot = qqplot_normal,
#        device = "jpg",
#        width = 10,
#        height = 5.6,
#        dpi = 600,
#        bg = "white")


## correlation between real and null pvalues
corr_table = bind_rows(null_pvals,
                       real_pvals) %>% 
  group_by(type, feature_name) %>% 
  arrange(pval) %>% 
  mutate(rank = row_number()) %>% 
  ungroup %>% 
  pivot_wider(names_from = type, values_from = pval)

# test difference with kolmogorov-smirnov
corr_table_ks = corr_table %>% 
  select(-rank) %>% 
  arrange(feature_name) %>% 
  split(as.factor(.$feature_name)) %>% 
  map(~dgof::ks.test(.x$`real features`, .x$`shuffled features`, alternative = "two.sided")) %>% 
  map(~broom::tidy(.x)) %>% 
  map(~pull(.x, p.value)) %>% 
  unlist %>% 
  data.frame %>% 
  rownames_to_column("feature_name") %>% 
  `colnames<-`(c("feature_name", "K-S p-value"))

write_tsv(corr_table_ks,
          "corr_table_ks.tsv")

# plot
corr_pvals_plot = ggplot(corr_table,
                         aes(x = `shuffled features`,
                             y = `real features`)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_text(data = corr_table_ks,
            aes(label = paste0("K-S p-value\n", round(`K-S p-value`, 2)), 
                x = 0.8, y = 0.2),
            size = 2) +
  facet_wrap(facets = vars(feature_name)) +
  ggtitle("p-values") +
  theme_classic() +
  theme(text = element_text(size = 6))
ggsave("7_corr_pvals_plot.jpg",
       plot = corr_pvals_plot,
       device = "jpg",
       width = 10,
       height = 5.6,
       dpi = 600,
       bg = "white")

# png("8_null_pvals_H3K36me3_roadmap.png")
# hist(null_pvals %>% 
#        filter(feature_name == 'H3K36me3_roadmap') %>% 
#        pull(pval), 
#      breaks = 250)
# dev.off()
# 
# png("9_real_pvals_H3K36me3_roadmap.png")
# hist(real_pvals %>% 
#        filter(feature_name == 'H3K36me3_roadmap') %>% 
#        pull(pval), 
#      breaks = 250)
# dev.off()
# 
# png("10_null_pvals_SUZ12.png")
# hist(null_pvals %>% 
#        filter(feature_name == 'SUZ12') %>% 
#        pull(pval), 
#      breaks = 250)
# dev.off()
# 
# png("11_real_pvals_SUZ12.png")
# hist(real_pvals %>% 
#        filter(feature_name == 'SUZ12') %>% 
#        pull(pval), 
#      breaks = 250)
# dev.off()
