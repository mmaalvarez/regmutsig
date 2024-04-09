library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")

regfeatures_list = read_tsv("../../../1st_half_genome/res/results_regional_feature_regressions_all_samples.tsv") %>% 
  filter(feature_type == "regional_feature") %>% 
  mutate(feature_name = gsub("_", " ", feature_name)) %>% 
  pull(feature_name) %>% 
  unique

res = read_tsv("rf_outputs/res.tsv") %>% 
  mutate(topredict = gsub(".*__resampled_", "Predictor resampled ", topredict),
         topredict = gsub("__|_", " ", topredict),
         topredict = gsub(" original", "", topredict),
         topredict = gsub("CI100%", "6 SD", topredict),
         topredict = gsub("2nd half ", "", topredict),
         topredict = factor(topredict,
                            levels = c(regfeatures_list,
                                       "SBS96",
                                       "Predictor 2nd genome half",
                                       "Predictor resampled 6 SD",
                                       paste0("Predictor resampled CI", c(99.99, 75, 50, 25, 5), "%"))),
         predictor = gsub("__", " - ", predictor),
         predictor = gsub("_", " ", predictor),
         predictor = gsub("1st half - | - original", "", predictor),
         predictor = gsub("regional feature", "Regional features", predictor)) %>% 
  arrange(predictor, topredict)


res_reg_feat_SBS96 = ggplot(res %>% 
                                  mutate(predictor_type = ifelse(predictor == "SBS96",
                                                                 "SBS96",
                                                                 "Regional feature"),
                                         topredict = gsub("1st genome half ", "", topredict)) %>%
                                  filter(!str_detect(topredict, "Predictor resampled")) %>%
                                  mutate(predictor = "1st genome half"),
                                aes(x = predictor, 
                                    y = topredict)) +
  geom_tile(aes(fill = R2)) +
  scale_fill_gradientn(colours = c('red','blue'),
                       limits = c(min(res$R2), max(res$R2))) +
  geom_label(aes(label = round(R2, 1)),
             size = 2,
             label.padding = unit(0.05, "lines"),
             label.size = 0) +
  facet_grid(cols = vars(predictor_type),
             switch = "x") +
  theme_classic() +
  xlab("Predictor (1st genome half)") +
  ylab("Predicted") +
  ggtitle("Inter-predictions of genome-split regression coefficients matrices") +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        strip.background = element_blank(),
        text = element_text(size = 10))
ggsave("res_reg_feat_SBS96.jpg",
       plot = res_reg_feat_SBS96,
       device = grDevices::jpeg,
       width = 5,
       height = 6,
       dpi = 600)

## scatterplot
res_reg_feat_SBS96_scatter = ggplot(res %>% 
                                      mutate(predictor_type = ifelse(predictor == "SBS96",
                                                                     "SBS96",
                                                                     "Regional feature"),
                                             topredict = gsub("1st genome half ", "", topredict)) %>%
                                      filter(!str_detect(topredict, "Predictor resampled")) %>%
                                      mutate(predictor = "1st genome half") %>% 
                                      pivot_wider(names_from = predictor_type,
                                                  values_from = R2) %>% 
                                      select(-predictor),
                                    aes(x = `Regional feature`, 
                                        y = SBS96)) +
  geom_point() +
  geom_abline(intercept = c(0,0), slope = 1, linetype = "dashed", color = "gray") +
  ggrepel::geom_text_repel(aes(label = topredict),
                           size = 2) +
  theme_bw() +
  xlab("R2 of a regional feature's ODD chromosomes SELF-predicting its EVEN chromosomes") +
  ylab("R2 of SBS96 ODD chromosomes predicting a regional feature's EVEN chromosomes") +
  theme(text = element_text(size = 10))
ggsave("res_reg_feat_SBS96_scatter.jpg",
       plot = res_reg_feat_SBS96_scatter,
       device = grDevices::jpeg,
       width = 15,
       height = 8.1,
       dpi = 600)
