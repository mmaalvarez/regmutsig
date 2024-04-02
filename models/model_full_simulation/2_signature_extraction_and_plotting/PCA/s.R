library(tidyverse)
library(data.table)
library(FactoMineR)
library(factoextra)
library(parallel)
library(cluster)
library(lsa)
library(cowplot)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("Position", "ggplot2")
conflict_prefer("cosine", "lsa")
conflict_prefer("clusterExport", "parallel")
conflict_prefer("clusterEvalQ", "parallel")
conflict_prefer("parLapply", "parallel")
conflict_prefer("parLapplyLB", "parallel")


results_regressions_basic_simulations = read_tsv("../../1_parser_and_regressions/res/results_regressions_basic_simulations.tsv") %>% 
  # not using conf. ints. for PCA
  filter(stat=="estimate") %>%
  select(-c(stat, theta, sample_type))
gc()

### Optional: IN CASE you want to extract specific parameters
# results_regressions_basic_simulations = results_regressions_basic_simulations %>%
#   filter(total_number_features == 8 &
#          number_tumors_per_feature == 25 &
#          total_number_normal_samples == 50 &
#          mean_fraction_of_nt_at_risk_mutated_in_target_bins == 2e-05 &
#          mean_fraction_of_nt_at_risk_mutated_in_nontarget_bins <= 1e-10)


## run PCA by group of parameters
pca_res = results_regressions_basic_simulations %>% 
  unite("parameters", total_number_features, number_tumors_per_feature, total_number_normal_samples, bin_size_bp, mean_fraction_of_nt_at_risk_mutated_in_target_bins, mean_fraction_of_nt_at_risk_mutated_in_nontarget_bins)

max_number_features = max(results_regressions_basic_simulations$total_number_features)

pca_res_split = pca_res %>% 
  group_split(parameters)
names(pca_res_split) = lapply(pca_res_split, function(x) unique(x$parameters))
pca_res_split = pca_res_split %>% 
  map(~pivot_wider(.x, names_from = factor, values_from = value)) %>% 
  map(~select(.x, -parameters)) %>%
  map(~column_to_rownames(.x, "sample_id")) %>% 
  map(~PCA(.x, graph = FALSE, ncp = max_number_features))

# create header for following loop
pca_res_final = pca_res_split[[length(names(pca_res_split))]]$ind$coord %>% 
  data.frame %>% 
  rownames_to_column("sample_id") %>% 
  #separate(sample_id, into = c("sample_id", "beta_10"), sep = "__") %>% 
  mutate(parameters = "") %>% 
  relocate(parameters, .after = sample_id) %>% 
  head(n=0)

for(parameters in names(pca_res_split)){
  pca_res_ind_coord = pca_res_split[[parameters]]$ind$coord %>% 
    data.frame %>% 
    mutate(parameters = parameters) %>% 
    rownames_to_column("sample_id") %>% 
    #separate(sample_id, into = c("sample_id", "beta_10"), sep = "__") %>% 
    relocate(parameters, .after = sample_id)
  pca_res_final = pca_res_final %>% 
    bind_rows(pca_res_ind_coord)
}

pca = pca_res_final %>% 
  separate(parameters, into = c('total_number_features', 'number_tumors_per_feature', 'total_number_normal_samples', 'bin_size_bp', 'mean_fraction_of_nt_at_risk_mutated_in_target_bins', 'mean_fraction_of_nt_at_risk_mutated_in_nontarget_bins'), sep = "_") %>% 
  rename_with(~str_replace(., 'Dim.', 'PC')) %>% 
  mutate(sample_id = gsub("_tumor.*", "-/-", sample_id),
         sample_id = gsub("normal.*", "Wild-type", sample_id)) %>% 
  rename("Deficient feature" = "sample_id",
         "Features" = "total_number_features",
         "Tumors×feature" = "number_tumors_per_feature",
         "Normals" = "total_number_normal_samples",
         "Bin size bp" = "bin_size_bp",
         "% nt@risk mut in target bins" = "mean_fraction_of_nt_at_risk_mutated_in_target_bins",
         "% mut rest" = "mean_fraction_of_nt_at_risk_mutated_in_nontarget_bins") %>% 
  mutate(`% nt@risk mut in target bins` = as.numeric(`% nt@risk mut in target bins`) * 100,
         `% mut rest` = as.numeric(`% mut rest`) * 100)



###########
### clustering to assess stability of clusters in each pca

# initialize cosine distance function which is used for clustering
cosineDist <- function(x) {
  mat = as.dist(t(1 - cosine(t(x))))
  return(mat)
}

# initialize function for clustering
nCPUs = 8
optimizeClustersPar = function(pca = NULL, maxNumClusters = max_number_features, threads = nCPUs) {
  
  cl = makeCluster(threads, useXDR = F) # 'makeCluster' refers to computer clusters!
  clusterEvalQ(cl, library(cluster)) # to be able to run "pam" ('cluster' refers to data clustering, not computer clusters!)  
  sigDiscoScores = data.table()
  
  pca_parameter_groups = pca %>% 
    unite("parameters", Features, `Tumors×feature`, Normals, `Bin size bp`, `% nt@risk mut in target bins`, `% mut rest`) %>%
    group_by(parameters, `Deficient feature`) %>% 
    mutate(sample_id = row_number()) %>% 
    ungroup %>% 
    unite("sample_id", `Deficient feature`, sample_id) %>% 
    dplyr::select(sample_id, parameters, PC1, PC2)
  
  for (parameter_group in unique(pca_parameter_groups$parameters)){
    
    pca_single_parameter_group = pca_parameter_groups %>% 
      filter(`parameters` == parameter_group) %>% 
      select(-parameters) %>% 
      column_to_rownames("sample_id")
    
    ## get matrix of cosine distances between samples' PC1 and PC2
    pca_parameter_group_Dists = cosineDist(pca_single_parameter_group)
    gc()
    
    k = str_extract(parameter_group, "[^_]+") %>% as.numeric
    
    # parallel code: in each thread, runs the "pam" clustering for a different # clusters, starting from the same matrix of cosine similarities
    # (note: Sanger (Alexandrov et al.) always run this with #clusters = #NMF factors, which may not be ideal)
    clusterExport(cl = cl, list("pca_parameter_group_Dists"), envir = environment()) # , envir=tand)
    pamOutputByIter = parLapplyLB(# load balancing
      cl = cl,
      X = k, # value of k in pam()
      fun = function(x) {
        set.seed(1)
        # partitioning (clustering) of the samples into k clusters "around medoids" (a more robust version of K-means), based on the cosine measures matrix
        pam(pca_parameter_group_Dists, k = x, diss = T)
      })
    gc()
    
    for (clusters in pamOutputByIter) {
      ## get clustering "stability": mean(clusters$silinfo$clus.avg.widths), i.e. average-across-K-clusters of the average silhouette widths per cluster
      # "silhouette width" of a sample (i.e. a given nFactor×run) i in a cluster is the "dissimilarity between i and the closest sample from an external cluster"
      sigDiscoScores = sigDiscoScores %>% 
        rbind(data.table(parameters = parameter_group, 
                         k = k,
                         avgClScore = mean(clusters$silinfo$clus.avg.widths), 
                         minClScore = min(clusters$silinfo$clus.avg.widths),
                         secondWorstClScore = clusters$silinfo$clus.avg.widths[order(clusters$silinfo$clus.avg.widths)][2],
                         avgPtScore = mean(silhouette(clusters)[, "sil_width"]), 
                         medianPtScore = median(silhouette(clusters)[, "sil_width"])))
      cat(sprintf("%02s\t%02d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                  parameter_group, k,
                  mean(clusters$silinfo$clus.avg.widths), 
                  min(clusters$silinfo$clus.avg.widths),
                  clusters$silinfo$clus.avg.widths[order(clusters$silinfo$clus.avg.widths)][2],
                  mean(silhouette(clusters)[, "sil_width"]), 
                  median(silhouette(clusters)[, "sil_width"])))
      gc()
    }
  }
  
  stopCluster(cl)
  #rm(maxNumClusters, threads, cl, pamOutputByIter, k )
  return(sigDiscoScores)
}

# run clustering
sigClusterScores = optimizeClustersPar(pca, 
                                       maxNumClusters = max_number_features) %>% 
  as_tibble %>% 
  dplyr::select(parameters, avgClScore) %>% 
  rename("Stability" = "avgClScore") %>% 
  separate(parameters, into = c('Features', 'Tumors×feature', 'Normals', 'Bin size bp', '% nt@risk mut in target bins', '% mut rest'), sep = "_")
  
write_tsv(sigClusterScores, "clusters_stabilities.tsv") # 'quality' of clustering 



####################################
### plot

jet.colors = colorRampPalette(c("gray", "red", "yellow", "green", "cyan", "blue", "magenta", "black"))

pca_grid_plot = ggplot(merge(pca, sigClusterScores), # add stability results
                       aes(x = PC1,
                           y = PC2)) +
  geom_point(aes(color = `Deficient feature`),
             size = 0.01,
             alpha = 0.8,
             shape = 16) +
  scale_color_manual(values = jet.colors(length(unique(pca$`Deficient feature`)))) +
  guides(color = guide_legend(override.aes = list(size=3, shape=16))) +
  ggh4x::facet_nested(rows = vars(`% nt@risk mut in target bins`, `% mut rest`),
                      cols = vars(`Features`, `Tumors×feature`, `Normals`),
                      scales = "free",
                      space = "free",
                      labeller = label_both) +
  geom_label(mapping = aes(label = round(Stability,1)),
             x = 0, y = -2.8,
             size = 1.5,
             label.padding = unit(units="lines", 0),
             label.size = unit(units="mm", 0)) +
  theme_classic() +
  theme(text = element_text(size = 3),
        axis.ticks = element_line(linewidth = 0.1),
        axis.line = element_line(linewidth = 0.1),
        strip.background = element_blank(),
        legend.title = element_blank())
ggsave("pca_grid.jpg",
       plot = pca_grid_plot,
       device = "jpg",
       width = 30,
       height = 11,
       dpi = 600)


### Optional: IN CASE you extracted specific parameters

#outliers_cutoff = 0.75

pc_x = "PC1"
pc_y = "PC2"

pca_top_sta_plot = ggplot(pca,
                          aes(x = !!sym(pc_x),
                              y = !!sym(pc_y))) +
  geom_point(aes(color = `Deficient feature`),
             size = 5,
             alpha = 1,
             shape = 16) +
  scale_color_manual(values = jet.colors(length(unique(pca$`Deficient feature`)))) +
  guides(color = guide_legend(override.aes = list(size=6, shape=16))) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  # ggrepel::geom_text_repel(data = pca[abs(pca[[pc_x]]) > outliers_cutoff*IQR(pca[[pc_x]]) | abs(pca[[pc_y]]) > outliers_cutoff*IQR(pca[[pc_y]]),],
  #                          aes(label = `Deficient feature`),
  #                          size = 3,
  #                          force = 20,
  #                          segment.size = 0.1,
  #                          min.segment.length = 0.001,
  #                          max.overlaps = 100000,
  #                          max.iter = 100000) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = "transparent"),
        text = element_text(size = 25))
ggsave(paste0("pca_top_stability_",pc_x,"_",pc_y,".jpg"),
       plot = pca_top_sta_plot,
       device = "jpg",
       width = 20,
       height = 10,
       dpi = 600)

vars = pca_res_split$`8_25_50_11718750_2e-05_1e-10`$var$coord %>%
  data.frame() %>%
  rename_with(~str_replace(., 'Dim.', 'PC')) %>% 
  rownames_to_column("vars") %>% 
  ggplot() +
  geom_segment(aes(x = 0, xend = !!sym(pc_x),
                   y = 0, yend = !!sym(pc_y)),
               arrow = arrow(length = unit(0.025,
                                           "npc"),
                             type = "open"),
               lwd = 0.5,
               linetype = "dashed") +
  ggrepel::geom_text_repel(aes(x = !!sym(pc_x),
                               y = !!sym(pc_y),
                               label = vars),
                           size = 3,
                           direction = "y",
                           vjust = 3,
                           force = 5,
                           segment.size = 0,
                           min.segment.length = 0,
                           max.overlaps = 1000000,
                           max.iter = 100000) +
  theme_nothing()
ggsave(paste0("vars_top_stability_", pc_x, "_", pc_y, ".jpg"),
       plot = vars,
       device = "jpg",
       width = 11,
       height = 5.6,
       dpi = 600,
       bg = "transparent")

scree = fviz_eig(pca_res_split$`8_25_50_11718750_2e-05_1e-10`,
                 ylim = c(0, round(max(pca_res_split$`8_25_50_11718750_2e-05_1e-10`$eig[,2]))),
                 geom = c("bar"),
                 addlabels = FALSE,
                 ncp = 8,
                 main = "",
                 ggtheme = theme_classic(base_size = 20),
                 xlab = "PC",
                 ylab = "% variance")
ggsave("scree_top_stability.jpg",
       plot = scree,
       device = "jpg",
       width = 10,
       height = 5.6,
       dpi = 600,
       bg = "transparent")
