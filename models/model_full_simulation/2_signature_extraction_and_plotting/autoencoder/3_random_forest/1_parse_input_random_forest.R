library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("Position", "ggplot2")


##### load and parse ae results

paths_encoded_layer_files = Sys.glob("../1_run_autoencoder/autoencoder_output/*/*.tsv")

encoded_layer_parameters = gsub("/encoded_layer.tsv", "", paths_encoded_layer_files) %>% 
  gsub(".*/", "", .) %>% 
  gsub("glorot_", "glorot.", .) %>% 
  gsub("cosine_", "cosine.", .) %>% 
  gsub("mean_squared_error", "mean.squared.error", .) %>% 
  data.frame %>% 
  separate(".", into = c('components', 'depth', 'hiddendim', 'batchsize', 'epochs', 'activation', 'weight_init', 'loss_function', 'learningrate', 'optimizer', 'epoch_min_val_loss', 'mean_sum_activity'), sep = "_") %>% 
  mutate_all(~gsub("-.*", "", .)) %>% 
  # keep track of order in paths_encoded_layer_files to then merge the encoded_layers components
  rownames_to_column("i") %>% 
  mutate(i = as.numeric(i)) %>% 
  relocate(i) %>% 
  select(-c(epoch_min_val_loss,mean_sum_activity))

encoded_layers = lapply(paths_encoded_layer_files,
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
  mutate(feature = factor(feature, levels = unique(.$feature), ordered = T)) %>% 
  arrange(feature) %>% 
  unite('parameters_group', components,depth,hiddendim,batchsize,epochs,activation,weight_init,loss_function,learningrate,optimizer, sep = "_") %>% 
  unite("sample_id", feature, Sample, sep = "_")


#################
## split encoded_layers table by i)parameters_group and ii)features, so each sample's simulated feature's status ('mut' (1) vs 'wt' (0)) is independently predicted within each autoencoder's parameters group

dir.create("rf_inputs")

for(single_parameters_group in unique(encoded_layers$parameters_group)){
  
  dir.create(paste0("rf_inputs/", single_parameters_group))
  
  encoded_layers_single_parameters_group = encoded_layers %>% 
    filter(parameters_group == single_parameters_group) %>% 
    select(-parameters_group)
  
  for(single_feature in paste0("feature", seq(1,8))){
    
    encoded_layers_single_parameters_group_single_feature = encoded_layers_single_parameters_group %>% 
      mutate(Phenotype = ifelse(str_detect(sample_id, paste0(single_feature,"_")), 1, 0),
             `AE component` = gsub("^", "ae", encoded_layer_component)) %>% 
      select(sample_id, Phenotype, `AE component`, Exposure) %>%
      pivot_wider(names_from = `AE component`, values_from = Exposure) %>% 
      # re-add the info on which samples were used for training and which for validation in the AE
      left_join(read_tsv("/g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model_full_simulation/2_signature_extraction_and_plotting/autoencoder/1_run_autoencoder/autoencoder_input/train_val_sets.tsv")) %>% 
      relocate(train_val) %>% 
      select(-sample_id)
    
    write_tsv(encoded_layers_single_parameters_group_single_feature,
              paste0("rf_inputs/", single_parameters_group, "/", single_feature, "_exposures.tsv"))
  }
}
