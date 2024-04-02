# Independent Component Analysis
# https://stats.stackexchange.com/questions/164872/the-independence-in-independent-component-analysis-intuitive-explanation
# https://www.statisticalaid.com/independent-component-analysis-ica-using-r/
# heuristic for finding if a feature gives clusters or not

library(tidyverse)
library(fastICA)
library(moments)


## load res
critical_value = 1.28 # for CI80%
# 1.96 for CI95%

coefficient_matrix = lapply(c("../../../1_parser_and_regressions/res/results_regional_feature_regressions_all_samples.tsv",
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
  filter(stat == 'estimate') %>% 
  select(Sample, feature_name, stat, val) %>%
  pivot_wider(names_from = feature_name, values_from = val) %>% 
  select(-stat) %>% 
  column_to_rownames("Sample") %>% 
  as.matrix
gc()


maxK = 12
res_list = list()
samples_vector = row.names(coefficient_matrix)
for(K in seq(2, maxK)){
  
  res = fastICA(X = coefficient_matrix, # rows samples, columns coefficients
                n.comp = K,
                alg.typ = "parallel",
                fun = "logcosh", 
                alpha = 1.0, 
                method = "C",
                row.norm = FALSE, 
                maxit = 200, 
                tol = 1e-04, 
                verbose = T)
  
  # calculate kurtosis of signatures (columns) in the output S ("estimated source") matrix
  kurt_vals = tibble("kurtosis" = kurtosis(res$S),
                     "Signature" = as.character(seq(1,K)))
  
  # plot
  S_parsed = res$S %>% 
    as_tibble %>% 
    bind_cols("Sample" = samples_vector) %>% 
    pivot_longer(cols = !matches("Sample"), names_to = "Signature", values_to = "Source (estimated)") %>% 
    mutate(K = paste0("K", K),
           Signature = gsub("V", "", Signature)) %>% 
    left_join(kurt_vals)
  
  res_list[[paste0("K: ", K)]] = S_parsed
}

res = bind_rows(res_list) %>% 
  mutate(Signature = factor(Signature, levels = seq(1,maxK)),
         K = factor(K, levels = paste0("K", seq(2,maxK)))) %>% 
  unite("K_Sig", K, Signature, sep = "\n", remove = F) %>% 
  mutate(K_Sig = factor(K_Sig, levels = unique(K_Sig)))

leptokurtic = res %>%
  # leptokurtic only
  select(Signature, K, K_Sig, kurtosis, `Source (estimated)`) %>%
  group_by(Signature, K, K_Sig, kurtosis) %>%
  slice_min(`Source (estimated)`) %>%
  distinct %>%
  filter(kurtosis > 3) %>%
  ungroup

total_leptokurtic_per_K = leptokurtic %>% 
  group_by(K) %>% 
  summarise(total_leptokurtic = n()) %>% 
  ungroup %>% 
  mutate(total_leptokurtic = paste0(total_leptokurtic, "/", gsub("K","",K), "=", round(total_leptokurtic/as.numeric(gsub("K","",K)),2), " leptokurtic,"))

total_possible_comparison_pairs = res %>% 
  # leptokurtic only
  filter(kurtosis>3) %>% 
  select(K_Sig, K) %>% 
  distinct %>% 
  group_by(K) %>% 
  summarise(total_possible_comparison_pairs = factorial(n()) / (2 * factorial(n() - 2))) %>% # n!/2×(n-2)!
  ungroup %>% 
  mutate(K_Sig = paste0(K,"\n", round(as.numeric(gsub("K", "", K))/2,0)))

mann_whitney_comparison_pairs = res %>% 
  # leptokurtic only
  filter(kurtosis>3) %>% 
  group_by(K) %>% 
  summarise(mw = list(pairwise.wilcox.test(`Source (estimated)`, K_Sig, p.adjust.method = "none", pool.sd = FALSE)$p.value %>% 
                        as.data.frame.table(responseName = "p.value"))) %>%
  unnest(mw) %>%
  filter(!is.na(p.value)) %>% 
  rename("K_Sig1" = "Var1",
         "K_Sig2" = "Var2") %>% 
  left_join(total_possible_comparison_pairs) %>% 
  left_join(total_leptokurtic_per_K) %>% 
  select(p.value, total_possible_comparison_pairs, K_Sig, total_leptokurtic) %>% 
  filter(p.value < 0.05) %>% 
  group_by(K_Sig, total_possible_comparison_pairs, total_leptokurtic) %>% 
  summarise(significant_comparison_pairs = n()) %>% 
  ungroup %>% 
  mutate("Fraction significant comparison pairs" = round(significant_comparison_pairs / total_possible_comparison_pairs, 2)) %>% 
  unite("Significant comparison pairs", significant_comparison_pairs, total_possible_comparison_pairs, sep = "/") %>% 
  unite("Significant comparison pairs", `Significant comparison pairs`, `Fraction significant comparison pairs`, sep = "=") %>% 
  mutate("Significant comparison pairs" = gsub("^", "of which signif. pairs\n", `Significant comparison pairs`),
         "Source (estimated)" = min(res$`Source (estimated)`)) %>% 
  unite("Significant comparison leptokurtic pairs", total_leptokurtic, `Significant comparison pairs`, sep = "\n") %>% 
  mutate(`Source (estimated)` = ifelse(str_detect(K_Sig, "K2"), `Source (estimated)` + 2,
                                       ifelse(str_detect(K_Sig, "K3"), `Source (estimated)` + 1, `Source (estimated)`)),
         K_Sig = gsub("K2\n1", "K2\n2", K_Sig))

# comparison_pairs = res %>% 
#   # leptokurtic only
#   filter(kurtosis>3) %>% 
#   select(K_Sig, K) %>% 
#   distinct %>% 
#   group_split(K) %>% 
#   map(~pull(.x, K_Sig)) %>% 
#   map(~combn(as.character(.), 2, simplify = F)) %>% 
#   flatten

jet_colors = colorRampPalette(c("red", "blue", "green", "magenta", "yellow", "cyan"))

S_plot = ggplot(res,
                aes(x = K_Sig,
                    y = `Source (estimated)`)) +
  coord_cartesian(ylim = c(min(res$`Source (estimated)`),
                           max(res$`Source (estimated)`))) +
  geom_violin(aes(fill = K)) +
  geom_boxplot(width = 0.1,
               outlier.size = 0.3,
               position=position_dodge(1)) +
  guides(fill = guide_legend(nrow=1, byrow=T)) +
  scale_fill_manual(values = jet_colors(length(levels(res$K)))) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "lightblue") +
  geom_label(data = leptokurtic,
             aes(label = round(kurtosis, 0)),
             size = 3,
             vjust = -2,
             label.padding = unit(0, "lines"),
             label.size = 0) +
  geom_text(data = mann_whitney_comparison_pairs,
            aes(label = `Significant comparison leptokurtic pairs`),
            size = 2) +
  # ggpubr::stat_compare_means(comparisons = comparison_pairs,
  #                            method = "wilcox.test",
  #                            method.args = list(alternative = "two.sided"),
  #                            hide.ns = T,
  #                            label = "p.signif",
  #                            step.increase = 0.001,
  #                            position="jitter",
  #                            vjust = 15,
  #                            bracket.size = 0,
  #                            symnum.args = list(cutpoints = c(0, 0.05, Inf), symbols = c(".", "ns"))) +
  xlab("Comparisons of K signatures (kurtosis label if >3; significance level '.' if Mann-Whitney p-value <0.05)") +
  ylab("Source (estimated) across samples") +
  theme_bw() +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(size = 5),
        legend.position = "top",
        legend.title = element_blank())
ggsave("S_plot.jpg",
       plot = S_plot,
       device = "jpg",
       width = 16,
       height = 9,
       dpi = 600)
