library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")

res = read_tsv("rf_outputs/res.tsv") %>% 
  #right_join(expand(., predictor, topredict)) %>% 
  mutate(topredict = gsub(".*__resampled_", "Predictor resampled ", topredict),
         topredict = gsub("__|_", " ", topredict),
         topredict = gsub("2nd half.*", "Predictor 2nd genome half", topredict),
         topredict = gsub(" original", "", topredict),
         topredict = gsub("1st half regional feature", "Regional features 1st genome half", topredict),
         topredict = gsub("CI100%", "6 SD", topredict),
         topredict = factor(topredict,
                            levels = c("Regional features 1st genome half",
                                       "Predictor 2nd genome half",
                                       "Predictor resampled 6 SD",
                                       paste0("Predictor resampled CI", c(99.99, 75, 50, 25, 5), "%"))),
         predictor = gsub("__", " - ", predictor),
         predictor = gsub("_", " ", predictor),
         predictor = gsub("1st half - | - original", "", predictor),
         predictor = gsub("regional feature", "Regional features", predictor)) %>% 
  arrange(predictor, topredict)

res_plot = ggplot(res,
                  aes(x = predictor, 
                      y = topredict)) +
  geom_tile(aes(fill = R2)) +
  scale_fill_gradientn(colours = c('red','blue')) +
  geom_label(aes(label = round(R2, 2)),
             size = 3) +
  facet_grid(cols = vars(predictor),
             scales = "free",
             space = "free") +
  theme_classic() +
  xlab("Predictor (1st genome half)") +
  ylab("Predicted") +
  ggtitle("Inter-predictions of genome-split regression coefficients matrices") +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text = element_blank(),
        strip.background = element_blank())
ggsave("res.jpg",
       plot = res_plot,
       device = grDevices::jpeg,
       width = 10.7,
       height = 6,
       dpi = 600)
