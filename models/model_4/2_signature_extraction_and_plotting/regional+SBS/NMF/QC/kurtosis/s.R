library(tidyverse)
library(moments)

## load res
critical_value = 1.28 # for CI80%
# 1.96 for CI95%

coefficient_matrix = lapply(c("../../../../../1_parser_and_regressions/res/results_regional_feature_regressions_all_samples.tsv",
                             "../../../../../../model_2/1_parser_and_regressions/res/results_SBS96_regressions_all_samples.tsv"),
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


range_nFact = seq(3, 12)
range_k = seq(3, 12)
optimal_nFact_k_list = expand.grid(range_nFact, # nFact
                                   range_k) %>% # K
  as.matrix %>% t %>% data.frame %>% as.list

res_list = list()
samples_vector = row.names(coefficient_matrix)

for(optimal_nFact_k in optimal_nFact_k_list){
  
  nFact = optimal_nFact_k[1]
  optimal_k = optimal_nFact_k[2] # the final number of signatures 
  
  if(!identical(Sys.glob(paste0("../../exposures_weights/fct", nFact, "_k", optimal_k, "_exposures.tsv")), character(0))){
  
    exposures = read_tsv(paste0("../../exposures_weights/fct", nFact, "_k", optimal_k, "_exposures.tsv"))
    
    # calculate kurtosis of samples exposures of NMF
    exposures_kurtosis = exposures %>% 
      select(Signature, Exposure) %>% 
      group_by(Signature) %>% 
      summarise(kurtosis = kurtosis(Exposure))
    
    # for plot
    exposures_kurtosis = exposures %>% 
      left_join(exposures_kurtosis) %>% 
      mutate(K = paste0("K=", optimal_k),
             nFact = paste0("nFact=", nFact))
    
    res_list[[paste0("nFact=", nFact, " K=", optimal_k)]] = exposures_kurtosis
  }
}

res = bind_rows(res_list) %>% 
  unite("nFact_K_Sig", nFact, K, Signature, sep = "\n", remove = F) %>% 
  mutate(nFact_K_Sig = factor(nFact_K_Sig, levels = unique(.$nFact_K_Sig)),
         K = factor(K, levels = paste0("K=", seq(max(range_k), min(range_k)))),
         nFact = factor(nFact, levels = unique(.$nFact)))

res_std_sig_names = res %>% 
  select(nFact, K, Signature) %>% 
  distinct %>% 
  group_by(nFact, K) %>% 
  mutate(Signature_std_name = paste0("Sig",row_number())) %>% 
  ungroup

res = left_join(res, 
                res_std_sig_names) %>% 
  mutate(Signature_std_name = factor(Signature_std_name,
                                     levels = paste0("Sig", seq(1, max(range_k)))))

leptokurtic = res %>%
  select(nFact, K, Signature_std_name, nFact_K_Sig, kurtosis, `Exposure`) %>%
  group_by(nFact, K, Signature_std_name, nFact_K_Sig, kurtosis) %>%
  slice_min(`Exposure`) %>%
  distinct %>%
  # leptokurtic only
  filter(kurtosis > 3) %>%
  ungroup

total_possible_comparison_pairs = res %>% 
  select(nFact_K_Sig, nFact, K) %>% 
  distinct %>% 
  group_by(nFact, K) %>% 
  summarise(total_possible_comparison_pairs = factorial(n()) / (2 * factorial(n() - 2))) %>% # n!/2×(n-2)!
  ungroup %>% 
  mutate(nFact_K_Sig = paste0(nFact, "\n", K,"\n", round(as.numeric(gsub("K=", "", K))/2,0)))

mann_whitney_comparison_pairs = res %>% 
  group_by(nFact, K) %>% 
  summarise(mw = list(pairwise.wilcox.test(`Exposure`, Signature_std_name, p.adjust.method = "none", pool.sd = FALSE)$p.value %>% 
                        as.data.frame.table(responseName = "p.value"))) %>%
  unnest(mw) %>%
  mutate(p.value = ifelse(is.na(p.value), 1, p.value)) %>% 
  left_join(total_possible_comparison_pairs) %>% 
  ungroup %>% 
  select(p.value, total_possible_comparison_pairs, nFact, K) %>% 
  filter(p.value < 0.05) %>% 
  group_by(nFact, K, total_possible_comparison_pairs) %>% 
  summarise(significant_comparison_pairs = n()) %>% 
  ungroup %>% 
  mutate(`Significant comparisons (%)` = (significant_comparison_pairs/total_possible_comparison_pairs)*100) %>% 
  select(nFact,K,`Significant comparisons (%)`)

res = left_join(res,
                mann_whitney_comparison_pairs)
leptokurtic = left_join(leptokurtic,
                        mann_whitney_comparison_pairs)

kurtosis_plot = ggplot(res,
                       aes(x = Exposure,
                           y = Signature_std_name,
                           fill = `Significant comparisons (%)`)) +
  coord_cartesian(xlim = c(0,
                           as.numeric(quantile(res$Exposure, 0.95)))) +
  geom_boxplot(outlier.alpha = 0,
               lwd=0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "lightblue") +
  geom_label(data = leptokurtic,
             aes(label = round(kurtosis, 0)),
             fill = "transparent",
             size = 2,
             nudge_x = as.numeric(quantile(res$Exposure, 0.95)),
             label.padding = unit(0, "lines"),
             label.size = 0) +
  geom_tile() +
  scale_fill_gradient2(low = "red", 
                       mid = "white",
                       high = "blue",
                       midpoint = median(mann_whitney_comparison_pairs$`Significant comparisons (%)`)) +
  facet_grid(cols = vars(nFact),
             rows = vars(K),
             space = "free", 
             scales = "free",
             switch = "both") +
  xlab("Exposure across samples (cap at overall 95% percentile)") +
  ylab("Signatures with kurtosis (>3) labels") +
  theme_bw() +
  theme(text = element_text(size = 10),
        legend.title = element_text(angle = 90))
ggsave("kurtosis_plot.jpg",
       plot = kurtosis_plot,
       device = "jpg",
       width = 16,
       height = 9,
       dpi = 600)
