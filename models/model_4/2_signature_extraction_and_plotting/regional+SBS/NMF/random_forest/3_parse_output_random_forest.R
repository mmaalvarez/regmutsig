library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")

res = read_tsv("rf_outputs/res.tsv") %>% 
  separate(NMF_parameters, into = c("nFact", "k"), sep = "_") %>%
  mutate(nFact = as.numeric(gsub("nFact-", "", nFact)),
         k = as.numeric(gsub("k-", "", k)),
         group = "Regional features + SBS") %>% 
  select(nFact, k, alteration, seed, oob, group) %>% 
  arrange(nFact, k, alteration, seed)

res_final_SBSonly = read_tsv("../../../../../model_2/2_signature_extraction_and_plotting/SBS/NMF/random_forest/rf_outputs/res_parsed.tsv")


### lots of parsing

res_final_tmp = res %>%
  bind_rows(res_final_SBSonly) %>% 
  mutate(group = factor(group, levels = c("Regional features + SBS", "SBS only"), ordered = T)) %>% 
  unite("k_nFact", k, nFact, sep = " N=") %>% 
  mutate(k_nFact = gsub("^", "K=", k_nFact)) 
# get highest OOB score per alteration and group
median_oob = res_final_tmp %>% 
  group_by(k_nFact, alteration, group) %>% 
  summarise(median_OOB = median(oob)) %>% 
  ungroup
# in case there is a tie between k_nFact combinations within the same alteration and group combination, keep the one with highest mean
mean_oob = res_final_tmp %>% 
  group_by(k_nFact, alteration, group) %>% 
  summarise(mean_OOB = mean(oob)) %>% 
  ungroup
  
res_final = merge(res_final_tmp, median_oob) %>% 
  merge(mean_oob) %>% 
  group_by(alteration, group) %>% 
  slice_max(median_OOB) %>% 
  slice_max(mean_OOB) %>% 
  ungroup %>% 
  mutate(alteration = factor(alteration))

# in case there is still a tie, remove the one with more negative skewness (i.e. tail towards negative oob)
skewness_oob = res_final %>% 
  select(alteration, group, k_nFact, oob) %>% 
  group_by(alteration, group, k_nFact) %>% 
  summarise(skewness = moments::skewness(oob)) %>% 
  # get the cases with >1 
  filter(n() >= 2)

# check whether there are also ties with SAME SKEWNESS...
same_skewness_oob = skewness_oob %>% 
  select(alteration, group, skewness) %>% 
  group_by(alteration, group, skewness) %>% 
  filter(n() >= 2) %>% 
  distinct %>% 
  mutate(same_skewness = "yes")

# if so, keep the one with the least number of clusters (K)
if(length(same_skewness_oob$alteration) >= 1){
  
  # store the ones to be removed in the final loop
  same_skewness_oob_to_rm = merge(skewness_oob,
                                  same_skewness_oob,
                                  all = T) %>% 
    mutate(same_skewness = ifelse(is.na(same_skewness), "no", same_skewness),
           K = as.numeric(gsub("K=| .*", "", k_nFact))) %>% 
    group_by(alteration, group, same_skewness) %>% 
    filter(same_skewness=="yes" & K != min(K)) %>% 
    group_by(alteration, group) %>% 
    select(alteration, group, k_nFact, skewness)
  
  skewness_oob = merge(skewness_oob,
                       same_skewness_oob,
                       all = T) %>% 
    mutate(same_skewness = ifelse(is.na(same_skewness), "no", same_skewness),
           K = as.numeric(gsub("K=| .*", "", k_nFact))) %>% 
    group_by(alteration, group, same_skewness) %>% 
    filter(same_skewness=="no" | K == min(K)) %>% 
    group_by(alteration, group) %>% 
    select(alteration, group, k_nFact, skewness)
} else {
  # just create a dummy one
  same_skewness_oob_to_rm = same_skewness_oob[FALSE, ]
}

# check whether there was a tie, loop until no ties left
while(length(skewness_oob$alteration) >= 2){
  
  # in case there is at least 1 tie, keep the more negatively skewed from each tie, to be removed
  more_neg_skewness_oob = skewness_oob %>% 
    slice_min(skewness) %>% 
    # append also the ones that have the same skewness and not the lowest K, to be removed as well
    bind_rows(same_skewness_oob_to_rm)
  
  # update skewness_oob without the current more neg. ones
  skewness_oob = skewness_oob %>% 
    filter(skewness != min(skewness)) %>% 
    # get the cases with still >1 
    group_by(alteration, group, k_nFact) %>% 
    filter(n() >= 2)
  
  res_final = res_final %>% 
    # append the skewness from the more negatively skewed ones (the ones that have the same skewness and not the lowest K) to table, to highlight them
    merge(more_neg_skewness_oob, all = T) %>% 
    # now keep only the ones without the skewness value
    filter(is.na(skewness)) %>% 
    select(-skewness)
}


## plotting

alterations_missing = res_final %>% 
  select(alteration, group) %>% 
  distinct() %>% 
  pull(alteration) %>% 
  table %>% 
  as.data.frame %>% 
  filter(Freq<=1) %>% 
  pull(".") %>% as.character

res_final_alterations = res_final %>% 
  filter(!alteration %in% alterations_missing) %>% 
  filter(alteration != 'UNG_GFP') %>%
  mutate(alteration = gsub("_KO|_and", "", alteration),
         alteration = gsub("_", "\n", alteration),
         alteration = factor(alteration))

# for labelling with Ks each box
alteration_group_KnFact_alterations = res_final_alterations %>% 
  select(alteration, group, k_nFact) %>% 
  distinct %>% 
  arrange(alteration, group) %>% 
  # dummy oob (for y position)
  mutate(oob = 0.7,
         # plot only the K
         K = gsub("K=| .*", "", k_nFact))

res_plot_alterations = ggplot(res_final_alterations,
                              aes(x = group, 
                                  y = oob)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  coord_cartesian(ylim = c(0.4, 1.03)) +
  geom_boxplot(aes(fill = group),
               size = 0.25,
               outlier.size = 0.1,
               fatten = 3) +
  scale_fill_manual(values = c("yellow", "blue")) +
  #geom_hline(yintercept = 0.9, linetype = "dashed", color = "darkgray", linewidth = 0.5) +
  ## mann whitney between models
  ggsignif::geom_signif(comparisons = combn(levels(res_final_alterations$group), 2, simplify = F),
                        map_signif_level = F, # asterisks if TRUE
                        test = "wilcox.test",
                        test.args = "two.sided") +  # or "less", "greater"
  geom_text(data = alteration_group_KnFact_alterations,
            aes(y=0.4,
                label = K),
            size = 3) +
  facet_grid(cols = vars(alteration)) +
  ggtitle("Predicted gene-KO, pathway alteration, specific treatment, or treatment type") +
  xlab("Number of signatures (K) that yields the highest median OOB score") +
  ylab("Out-of-bag score (cap at 0.4)") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.size = unit(1.5, 'cm'),
        legend.text = element_text(size = 20),
        text = element_text(size = 20),
        axis.text.x = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5),
        strip.text = element_text(size = 8, angle = 90),
        strip.background.x = element_blank())
ggsave("res.jpg",
       plot = res_plot_alterations,
       device = grDevices::jpeg,
       width = 16,
       height = 9,
       dpi = 600)


############################################

### explore some interesting alterations
# change the "alteration" and "K\\=" parts

res_final_alterations_interesting = res_final_alterations %>% 
  filter(
         (alteration == "Heterocyclic\namines" & ((group == "Regional features + SBS" & str_detect(k_nFact, "K\\=8"))
                                                  | 
                                                  (group == "SBS only" & str_detect(k_nFact, "K\\=11"))
                                                  )
          )
         |
         (alteration == "TLS" & ((group == "Regional features + SBS" & str_detect(k_nFact, "K\\=11"))
                                 | 
                                 (group == "SBS only" & str_detect(k_nFact, "K\\=6"))
                                 )
          )
         ) %>% 
  rename("altered_pathway_or_treatment_type" = "alteration") %>% 
  select(group, k_nFact, altered_pathway_or_treatment_type, mean_OOB, median_OOB) %>% 
  distinct %>% 
  arrange(group, k_nFact, altered_pathway_or_treatment_type)


#### according to res_final_alterations_interesting, extract the 2 winning combinations of nFact-k (1 for reg+SBS and another for only SBS) for each of the winning alterations
exposures_table1 = read_tsv(paste0("../exposures_weights/fct6_k8_exposures.tsv")) %>% 
  filter(altered_pathway_or_treatment_type == "Heterocyclic amines") %>% 
  mutate(group = "Regional features + SBS",
         k_nFact = "K=8 N=6",
         Signature = gsub(" - ", " ", Signature))
exposures_table2 = read_tsv(paste0("../../../SBS/NMF/exposures_weights/fct6_k11_exposures.tsv")) %>% 
  filter(altered_pathway_or_treatment_type == "Heterocyclic amines") %>% 
  mutate(group = "SBS only",
         k_nFact = "K=11 N=6",
         Signature = gsub(" - ", " ", Signature))
exposures_table3 = read_tsv(paste0("../exposures_weights/fct11_k11_exposures.tsv")) %>% 
  filter(altered_pathway_or_treatment_type == "TLS") %>% 
  mutate(group = "Regional features + SBS",
         k_nFact = "K=11 N=11",
         Signature = gsub(" - ", " ", Signature))
exposures_table4 = read_tsv(paste0("../../../SBS/NMF/exposures_weights/fct12_k6_exposures.tsv")) %>% 
  filter(altered_pathway_or_treatment_type == "TLS") %>% 
  mutate(group = "SBS only",
         k_nFact = "K=6 N=12",
         Signature = gsub(" - ", " ", Signature))
# bindrow them
exposures_table_interesting = bind_rows(exposures_table1,
                                        exposures_table2,
                                        exposures_table3,
                                        exposures_table4) %>% 
  select(group, k_nFact, altered_pathway_or_treatment_type, alteration, Sample, dataset, Signature, Exposure) %>% 
  distinct %>% 
  arrange(group, k_nFact, altered_pathway_or_treatment_type, alteration, dataset, Sample, Signature)


##### add to each signature the weights from each regional feature and SBS
weights_table1 = read_tsv(paste0("../exposures_weights/fct6_k8_weights.tsv")) %>% 
  select(Signature, Stability, feature, feature_group, Weight, `Max. sim. COSMIC`) %>% 
  mutate(group = "Regional features + SBS",
         k_nFact = "K=8 N=6")
weights_table2 = read_tsv(paste0("../../../SBS/NMF/exposures_weights/fct6_k11_weights.tsv")) %>% 
  rename("feature" = "SBS96") %>% 
  mutate(feature_group = "SBS") %>% 
  select(Signature, Stability, feature, feature_group, Weight, `Max. sim. COSMIC`) %>% 
  mutate(Signature = gsub(" - ", " ", Signature),
         group = "SBS only",
         k_nFact = "K=11 N=6",
         Stability = as.numeric(gsub("Stability: ",
                                     "",
                                     Stability)))
weights_table3 = read_tsv(paste0("../exposures_weights/fct11_k11_weights.tsv")) %>% 
  select(Signature, Stability, feature, feature_group, Weight, `Max. sim. COSMIC`) %>% 
  mutate(group = "Regional features + SBS",
         k_nFact = "K=11 N=11")
weights_table4 = read_tsv(paste0("../../../SBS/NMF/exposures_weights/fct12_k6_weights.tsv")) %>% 
  rename("feature" = "SBS96") %>% 
  mutate(feature_group = "SBS") %>% 
  select(Signature, Stability, feature, feature_group, Weight, `Max. sim. COSMIC`) %>% 
  mutate(Signature = gsub(" - ", " ", Signature),
         group = "SBS only",
         k_nFact = "K=6 N=12",
         Stability = as.numeric(gsub("Stability: ",
                                     "",
                                     Stability)))
# bindrow them
weights_table_interesting = bind_rows(weights_table1,
                                      weights_table2,
                                      weights_table3,
                                      weights_table4) %>% 
  select(group, k_nFact, feature, feature_group, Signature, Weight, Stability, `Max. sim. COSMIC`) %>% 
  distinct %>% 
  arrange(group, k_nFact, feature_group, feature, Signature) %>% 
  select(-c(feature_group)) %>% 
  # spread names as columns
  pivot_wider(names_from = feature, values_from = Weight)


### finish main full table and print it out
exposures_weights_interesting = exposures_table_interesting %>% 
  left_join(weights_table_interesting) 

write_tsv(exposures_weights_interesting,
          "exposures_weights_interesting.tsv")


### keep exploring the most important columns (exposures and feature weights)
exposure_thr = 5
weight_thr = 1e-2 # make sure it's lower than exposure_thr
filtered_table = exposures_weights_interesting %>% 
  # remove SBS only rows
  filter(group != "SBS only") %>% 
  # remove SBS columns, and other info
  select(-c(group, k_nFact, dataset, Stability, `Max. sim. COSMIC`, contains(">"))) %>% 
  # keep only signature×sample rows in which the exposure of that sample to that signature is above some threshold
  filter(Exposure >= exposure_thr) %>% 
  # keep only column weights in which the abs(weight) is larger than a threshold
  select_if(~ is.character(.) || (is.numeric(.) && any(abs(.) > weight_thr)))
