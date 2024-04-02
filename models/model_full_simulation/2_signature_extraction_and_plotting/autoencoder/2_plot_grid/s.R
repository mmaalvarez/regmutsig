library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("Position", "ggplot2")


##### load and parse autoencoder results

paths_encoded_layer_files = Sys.glob("../1_run_autoencoder/autoencoder_output/*/*.tsv")

encoded_layer_parameters = gsub("/encoded_layer.tsv", "", paths_encoded_layer_files) %>% 
  gsub(".*/", "", .) %>% 
  gsub("glorot_", "glorot.", .) %>% 
  gsub("cosine_", "cosine.", .) %>% 
  gsub("mean_squared_error", "mean.squared.error", .) %>% 
  data.frame %>% 
  separate(".", into = c('validation', 'components', 'depth', 'hiddendim', 'batchsize', 'epochs', 'activation', 'weight_init', 'loss_function', 'learningrate', 'optimizer', 'mean_sum_activity', 'seed'), sep = "_") %>% 
  mutate_all(~gsub("-.*", "", .)) %>% 
  # keep track of order in paths_encoded_layer_files to then merge the encoded_layers components
  rownames_to_column("i") %>% 
  mutate(i = as.numeric(i)) %>% 
  relocate(i)
  
encoded_layers_avg = lapply(paths_encoded_layer_files,
                        read_tsv) %>% 
  map2(seq_along(.), ~mutate(.x, i = .y)) %>% 
  bind_rows %>% 
  left_join(encoded_layer_parameters) %>% 
  select(-i) %>% 
  pivot_longer(cols = starts_with("ae"), names_to = "encoded_layer_component", values_to = "Exposure") %>% 
  mutate(encoded_layer_component = gsub("ae", "", encoded_layer_component),
         encoded_layer_component = factor(encoded_layer_component, levels = as.character(seq(1,length(unique(.$encoded_layer_component)))), ordered = T),
         Sample = factor(Sample, levels = sort(unique(.$Sample))),
         # make that hiddendim be the multiple of components that was used to define the hidden layers, rather than the resulting n neurons thereof
         hiddendim = as.numeric(hiddendim)/as.numeric(components)) %>% 
  # some only have 7 comp, others only 8
  filter(!is.na(`Exposure`)) %>% 
  # have a "feature" column
  separate("Sample", into = c("feature", "Sample"), sep = "_") %>% 
  mutate(Sample = ifelse(is.na(Sample),
                         feature,
                         Sample),
         feature = gsub("normal.*", "normal", feature)) %>% 
  # average exposures for all samples that have the same feature "mutated", or are all normal, and that have the same parameters 
  group_by(feature, validation, components,depth, hiddendim , batchsize , epochs ,activation ,weight_init, loss_function,learningrate,optimizer, seed, encoded_layer_component) %>% 
  summarize(`Mean exposure` = mean(Exposure)) %>% 
  ungroup

write_tsv(encoded_layers_avg, "encoded_layers_avg.tsv")


### plotting
encoded_layers_plot = ggplot(encoded_layers_avg, 
                             aes(x = feature,
                                 y = encoded_layer_component)) +
  geom_tile(aes(fill = `Mean exposure`)) +
  scale_fill_gradientn(colours = c('white','red')) + 
  ggh4x::facet_nested(rows = vars(depth,components,batchsize,epochs,learningrate),
                      # for depth==1, hiddendim acts as replicates (diff seed)
                      cols = vars(activation,optimizer,loss_function,validation,weight_init,hiddendim),
                      scales = "free",
                      space = "free",
                      labeller = label_both) +
  xlab("Pooled samples for each feature") + 
  ylab("Encoded layer component") +
  theme_bw() +
  theme(text = element_text(size = 7),
        axis.ticks = element_line(linewidth = 0.1),
        axis.line = element_line(linewidth = 0.1),
        axis.text.x = element_text(size = 2, angle=45, hjust=1),
        axis.text.y = element_text(size = 2),
        legend.position = "top")
ggsave(paste0("encoded_layers_simulations_grid.jpeg"),
       plot = encoded_layers_plot,
       device = "jpeg",
       width = 25,
       height = 14,
       dpi = 600)


####################################
### plot specific parameters

## fix clearly better parameters
encoded_layers_fixed_relu_adam_mse = encoded_layers_avg %>% 
  filter(activation == 'relu' &
         optimizer == 'Adam' &
         loss_function == 'mean.squared.error')

encoded_layers_fixed_relu_adam_mse_plot = ggplot(encoded_layers_fixed_relu_adam_mse, 
                                        aes(x = feature,
                                            y = encoded_layer_component)) +
  geom_tile(aes(fill = `Mean exposure`)) +
  scale_fill_gradientn(colours = c('white','red')) + 
  ggh4x::facet_nested(rows = vars(depth,batchsize,epochs,learningrate),
                      # for depth==1, hiddendim acts as replicates (diff seed)
                      cols = vars(components,validation,weight_init,hiddendim),
                      scales = "free",
                      space = "free",
                      labeller = label_both) +
  xlab("Pooled samples for each feature") + 
  ylab("Encoded layer component") +
  theme_bw() +
  theme(text = element_text(size = 10),
        axis.ticks = element_line(linewidth = 0.5),
        axis.line = element_line(linewidth = 0.5),
        axis.text.x = element_text(size = 4, angle=45, hjust=1),
        axis.text.y = element_text(size = 4),
        legend.position = "top")
ggsave(paste0("encoded_layers_fixed_relu_adam_mse.jpeg"),
       plot = encoded_layers_fixed_relu_adam_mse,
       device = "jpeg",
       width = 24.2,
       height = 14,
       dpi = 600)


## fix more parameters
encoded_layers_fixed_relu_adam_mse_depth2_k8_9 = encoded_layers_avg %>% 
  filter(activation == 'relu' &
           optimizer == 'Adam' &
           loss_function == 'mean.squared.error' &
           depth == '2' &
           components %in% c('8','9'))

encoded_layers_fixed_relu_adam_mse_depth2_k8_9_plot = ggplot(encoded_layers_fixed_relu_adam_mse_depth2_k8_9, 
                                                             aes(x = feature,
                                                                 y = encoded_layer_component)) +
  geom_tile(aes(fill = `Mean exposure`)) +
  scale_fill_gradientn(colours = c('white','red')) + 
  ggh4x::facet_nested(rows = vars(components,batchsize,epochs,learningrate),
                      # for depth==1, hiddendim acts as replicates (diff seed)
                      cols = vars(validation,weight_init,hiddendim),
                      scales = "free",
                      space = "free",
                      labeller = label_both) +
  xlab("Pooled samples for each feature") + 
  ylab("Encoded layer component") +
  theme_bw() +
  theme(text = element_text(size = 10),
        axis.ticks = element_line(linewidth = 0.5),
        axis.line = element_line(linewidth = 0.5),
        axis.text.x = element_text(size = 4, angle=45, hjust=1),
        axis.text.y = element_text(size = 4),
        legend.position = "top")
ggsave(paste0("encoded_layers_fixed_relu_adam_mse_depth2_k8_9.jpeg"),
       plot = encoded_layers_fixed_relu_adam_mse_depth2_k8_9_plot,
       device = "jpeg",
       width = 24.2,
       height = 14,
       dpi = 600)


## fix more parameters
encoded_layers_fixed_relu_adam_mse_depth2_k8_9_hidden1 = encoded_layers_avg %>% 
  filter(activation == 'relu' &
           optimizer == 'Adam' &
           loss_function == 'mean.squared.error' &
           depth == '2' &
           components %in% c('8','9') &
           hiddendim == '1')

encoded_layers_fixed_relu_adam_mse_depth2_k8_9_hidden1_plot = ggplot(encoded_layers_fixed_relu_adam_mse_depth2_k8_9_hidden1, 
                                                             aes(x = feature,
                                                                 y = encoded_layer_component)) +
  geom_tile(aes(fill = `Mean exposure`)) +
  scale_fill_gradientn(colours = c('white','red')) + 
  ggh4x::facet_nested(rows = vars(components,batchsize,epochs),
                      # for depth==1, hiddendim acts as replicates (diff seed)
                      cols = vars(validation,weight_init,learningrate),
                      scales = "free",
                      space = "free",
                      labeller = label_both) +
  xlab("Pooled samples for each feature") + 
  ylab("Encoded layer component") +
  theme_bw() +
  theme(text = element_text(size = 10),
        axis.ticks = element_line(linewidth = 0.5),
        axis.line = element_line(linewidth = 0.5),
        axis.text.x = element_text(size = 4, angle=45, hjust=1),
        axis.text.y = element_text(size = 4),
        legend.position = "top")
ggsave(paste0("encoded_layers_fixed_relu_adam_mse_depth2_k8_9_hidden1.jpeg"),
       plot = encoded_layers_fixed_relu_adam_mse_depth2_k8_9_hidden1_plot,
       device = "jpeg",
       width = 24.2,
       height = 14,
       dpi = 600)