library(tidyverse)

res = read_tsv("res.tsv") %>% 
  select(-contains("_importance")) %>% 
  pivot_longer(cols = -c(autoencoder_parameters, feature, test_size, seed),
               names_to = "metric",
               values_to = "value") %>% 
  mutate(metric = factor(metric, levels = rev(unique(.$metric)))) %>% 
  arrange(autoencoder_parameters, feature, test_size, seed) %>% 
  separate(autoencoder_parameters, into = c("components","depth","hiddendim","batchsize","epochs","activation","weight_init","loss_function","learningrate","optimizer"), sep = "_") %>% 
  select(-c(activation, weight_init, loss_function, optimizer))

res_plot = ggplot(res,
                  aes(x = metric,
                      y = value)) +
  geom_boxplot() +
  ggh4x::facet_nested(rows = vars(depth,hiddendim,batchsize,epochs,learningrate), 
                      cols = vars(components,feature,test_size),
                      scales = "free", space = "free",
                      labeller=label_both) +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "lightblue") +
  xlab("") +
  ylab("Metric value") +
  ggtitle(paste0("Results for 100 random forests: prediction of DNA repair phenotype (wt or deficient) using the exposures of the AE encoded components")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, vjust = 1))
ggsave("res.jpg",
       plot = res_plot,
       device = "jpeg",
       width = 27,
       height = 15,
       dpi = 600)
