library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")

regfeatures_list = read_tsv("../../1st_half_genome/res/results_regional_feature_regressions_all_samples.tsv") %>% 
  filter(feature_type == "regional_feature") %>% 
  mutate(feature_name = gsub("_", " ", feature_name)) %>% 
  pull(feature_name) %>% 
  unique

res = read_tsv("rf_outputs/res.tsv") %>% 
  mutate(topredict = gsub(".*__resampled_", "Predictor resampled ", topredict),
         topredict = gsub("__|_", " ", topredict),
         topredict = gsub("2nd half.*", "Predictor 2nd genome half", topredict),
         topredict = gsub(" original", "", topredict),
         topredict = gsub("1st half ", "1st genome half ", topredict),
         topredict = gsub("CI100%", "6 SD", topredict),
         topredict = factor(topredict,
                            levels = c(paste0("1st genome half ", regfeatures_list),
                                       "Predictor 2nd genome half",
                                       "Predictor resampled 6 SD",
                                       paste0("Predictor resampled CI", c(99.99, 75, 50, 25, 5), "%"))),
         predictor = gsub("__", " - ", predictor),
         predictor = gsub("_", " ", predictor),
         predictor = gsub("1st half - | - original", "", predictor),
         predictor = gsub("regional feature", "Regional features", predictor)) %>% 
  arrange(predictor, topredict)

res_reg_feat_predictor = ggplot(res %>% 
                                  filter(predictor != "SBS96"),
                                aes(x = predictor, 
                                    y = topredict)) +
  geom_tile(aes(fill = R2)) +
  scale_fill_gradientn(colours = c('red','blue'),
                       limits = c(min(res$R2), max(res$R2))) +
  geom_label(aes(label = round(R2, 1)),
             size = 1.1,
             label.padding = unit(0.05, "lines"),
             label.size = 0) +
  facet_grid(cols = vars(predictor),
             scales = "free",
             space = "free") +
  theme_classic() +
  xlab("Predictor (1st genome half)") +
  ylab("Predicted") +
  ggtitle("Inter-predictions of genome-split regression coefficients matrices") +
  theme(plot.title = element_text(hjust = 1),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(hjust = 1),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")
ggsave("res_reg_feat_predictor.jpg",
       plot = res_reg_feat_predictor,
       device = grDevices::jpeg,
       width = 10.7,
       height = 6,
       dpi = 600)

res_SBS96_predictor = ggplot(res %>% 
                              filter(predictor == "SBS96"),
                      aes(x = predictor, 
                          y = topredict)) +
  geom_tile(aes(fill = R2)) +
  scale_fill_gradientn(colours = c('red','blue'),
                       limits = c(min(res$R2), max(res$R2))) +
  geom_label(aes(label = round(R2, 1)),
             size = 1.1,
             label.padding = unit(0.1, "lines"),
             label.size = 0) +
  facet_grid(cols = vars(predictor),
             scales = "free",
             space = "free") +
  theme_classic() +
  xlab("") +
  ylab("") +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_blank(),
        strip.background = element_blank())
ggsave("res_SBS96_predictor.jpg",
       plot = res_SBS96_predictor,
       device = grDevices::jpeg,
       width = 4,
       height = 6,
       dpi = 600)
