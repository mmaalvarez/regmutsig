library(tidyverse)
library(data.table)
library(parallel)
library(msm) # rtnorm
#install.packages('terra', repos='https://rspatial.r-universe.dev')
library(NMF) # for NMF to determine optimal k signatures
#ENMTools::install.extras("NMF")
library(lsa)
library(cowplot)
#library(cluster)
library(ggrepel)
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
conflict_prefer("nmf", "NMF")
conflict_prefer("Position", "ggplot2")



## results regressions

coefficient_table = read_tsv("../../1_parser_and_regressions/res/results_regressions_basic_simulations.tsv") %>%
  ## keep parameters of the simulated data whose regression results you want to run NMF on
  filter(total_number_features == 8 & 
         number_tumors_per_feature == 20 &
         total_number_normal_samples == 100 &
         mean_fraction_of_nt_at_risk_mutated_in_target_bins == 2e-5 &
         mean_fraction_of_nt_at_risk_mutated_in_nontarget_bins == 5e-8) %>% 
  select(sample_id, factor, stat, value) %>% 
  pivot_wider(names_from = stat, values_from = value) %>% 
  ## DONT DO THIS (removes signal): at samples with any estimate(s) >~|10|, convert their -culprit- estimates and CIs to 0 for the NMF, as these would mean that the regression died
  # rowwise %>% 
  # mutate(beta10 = ifelse(abs(estimate)>=10, "yes", "no"),
  #        across(starts_with("estimate") | starts_with("conf."),
  #               ~ifelse(beta10 == "yes", 0, .))) %>% 
  arrange(sample_id, factor) %>% 
  group_by(sample_id, factor)


#### permutations

## Parameters and initializing of some objects
totalNumIters = 100 #1000
coefficient_Resamp = list()
set.seed(1)

## Generate matrices resampling betas from their CI95% distributions (instead of UPmultinomial)
resample_from_CI = function(coefficient_table){
  summarise(coefficient_table, resampled_estimate = rtnorm(n = 1,
                                                           mean = estimate,
                                                           sd = 1, 
                                                           lower = conf.low,
                                                           upper = conf.high))}
for (nIter in 1:totalNumIters) {
  print(paste0("nIter ", nIter))
  # for each sample (row) resample coefficients from CI95% distrs.
  coefficient_Resamp[[nIter]] = resample_from_CI(coefficient_table) %>% 
    mutate("nIter" = nIter)
  gc()
}

# bind all permuted tables
coefficient_Resamp = bind_rows(coefficient_Resamp)

dir.create("heatmap_Ks")
write_tsv(coefficient_Resamp, paste0("heatmap_Ks/permuted_coefficients_", totalNumIters, "iters.tsv"))
#coefficient_Resamp = read_tsv(paste0("heatmap_Ks/permuted_coefficients_", totalNumIters, "iters.tsv"))



#####################################################################
#### NMF
#####################################################################


#### evaluate best combination of n variables and k signatures

coefficient_Resamp = coefficient_Resamp %>% 
  ungroup %>% 
  pivot_wider(names_from = factor, values_from = resampled_estimate) %>% 
  unite("sample_id", sample_id, nIter, sep = "__")

# a) keep positive coefficients, and convert negative to zero 
coefficient_Resamp_posmatrix = coefficient_Resamp %>% 
  mutate_if(is.numeric,
            ~if_else(.<=0, 0, .))
# b) convert positive coefficients to zero, and convert negative to positive
coefficient_Resamp_negmatrix = coefficient_Resamp %>% 
  mutate_if(is.numeric,
            ~if_else(.>=0, 0, abs(.)))
# merge converted coefficients into the NMF input
coefficient_Resamp = merge(coefficient_Resamp_posmatrix,
                           coefficient_Resamp_negmatrix,
                           by = "sample_id",
                           suffixes = c("_poscoeff", "_negcoeff")) %>%
  mutate(sample_id = gsub("__muts", "_muts", sample_id)) %>% 
  separate(sample_id, into = c("id", "nIter"), sep = "__", remove = F) %>% 
  mutate(sample_id = gsub("__", "_nIter", sample_id)) %>% 
  column_to_rownames("sample_id") %>% 
  select(-id) %>% 
  split(., f = .$nIter) %>% 
  map(.f = list(. %>% data.matrix))



## Run NMF for each resampled coefficients matrix generated in the previous step

## Parameters and initializing of some objects
set.seed(1)
maxK = 10
# resampled NMF input matrices, generated just once and then re-used for any possible # factors/clusters
nmfHmatAllByFact = list()
nmfWmatAllByFact = list()
nCPUs = 8

## "cl" is only used if using parLapply (creates a set of copies of R running in parallel and communicating over sockets, slower?) instead of mclapply (parallel processing):
#cl = makeCluster(nCPUs)
#clusterExport(cl, list("coefficient_Resamp"))
#clusterEvalQ(cl, library(NMF))

## extract nFact factors...
for(nFact in 2:maxK){
  
  cat(sprintf("Running NMF: nFact %d (all iters)\n", nFact))
  
  # ...from each of the permuted matrices
  nmfOutputByIter = mclapply(coefficient_Resamp,
                             mc.cores = nCPUs,
                             function(x, nFact){
                               # store colnames
                               orig_colnames_order = colnames(x)[-1] %>% 
                                 gsub("\\&|\\-", "\\.", .)
                               # NAs to 0
                               x[is.na(x)] <- 0
                               # store nIter
                               nIter = x[1,1]
                               x = x[,-1]
                               ## cannot handle all-zero columns (such as MSH6-neg, since all coeff. are +...)
                               # store their colnames
                               all_zero_columns = data.frame(x) %>% 
                                                    select_if(colSums(.) <= 0) %>% 
                                                    colnames
                               # remove them from matrix before NMF
                               x_nozerocols = x[, colSums(x) > 0]
                               ## cannot handle all-zero rows either, remove them from matrix before NMF (not re-added afterwards)
                               x_nozerocolsrows = x_nozerocols[rowSums(x_nozerocols) > 0,]
                               ### NMF
                               nmf_run = NMF::nmf(x_nozerocolsrows, 
                                                  rank = nFact, 
                                                  maxIter = 10000, 
                                                  seed = nIter)
                               # re-add the all-zero columns as zeros
                               all_zero_columns_matrix = data.frame(matrix(rep(0, len = length(all_zero_columns)), 
                                                                           ncol = length(all_zero_columns), 
                                                                           nrow = nrow(nmf_run@fit@H))) %>% 
                                 `colnames<-`(all_zero_columns) %>% 
                                 as.matrix
                               full_matrix = cbind(nmf_run@fit@H, all_zero_columns_matrix) %>% 
                                 data.frame
                               ## when making it data.frame the '(' '>' ')' '×' symbols become '.', so revert this
                               #full_matrix_colnames = colnames(full_matrix) %>% 
                               # str_replace(., 
                               #             "([A,C,G,T])\\.([C,T])\\.([A,C,G,T])\\.([A,C,G,T])_", 
                               #             "\\1(\\2>\\3)\\4_") %>% 
                               # str_replace(., "\\.", "\\×")
                               #colnames(full_matrix) = full_matrix_colnames
                               full_matrix = full_matrix %>% 
                                 # reorder column names to be the same order as before
                                 select(all_of(orig_colnames_order)) %>%
                                 as.matrix
                               nmf_run@fit@H = full_matrix
                               nmf_run},
                             nFact)
  idString = sprintf("nFact=%03d", nFact)
  nmfHmatAllByFact[[idString]] = matrix(nrow = 0, ncol = ncol(nmfOutputByIter[[1]])); # columns = original features...
  nmfWmatAllByFact[[idString]] = matrix(nrow = 0, ncol = nrow(nmfOutputByIter[[1]])); # columns = samples...
  for (nIter in 1:totalNumIters) {
    nmfOutput = nmfOutputByIter[[nIter]]
    # weights (H)
    rownames(nmfOutput@fit@H) = sprintf("run%03d_fact%02d", nIter, 1:nFact)
    # exposures ("basis", W)
    colnames(nmfOutput@fit@W) = sprintf("run%03d_fact%02d", nIter, 1:nFact)
    nmfHmatAllByFact[[idString]] = nmfHmatAllByFact[[idString]] %>% rbind(nmfOutput@fit@H) # same thing that you get using the coef(NMFfit) command in R   -> loadings
    nmfWmatAllByFact[[idString]] = nmfWmatAllByFact[[idString]] %>% rbind(t(nmfOutput@fit@W)) # same thing that you get using the basis(NMFfit) command in R  -> transformed data
    gc()
  }
  gc()
}
#stopCluster(cl)
rm(idString, nmfOutput, nFact)

## Save clustering results
saveRDS(nmfHmatAllByFact, paste0("heatmap_Ks/", "Hmatrix.Rda")) # signatures
#nmfHmatAllByFact = readRDS(paste0("heatmap_Ks/", "Hmatrix.Rda"))
saveRDS(nmfWmatAllByFact, paste0("heatmap_Ks/", "Wmatrix.Rda")) # exposures
#nmfWmatAllByFact = readRDS(paste0("heatmap_Ks/", "Wmatrix.Rda"))


## abs([results of pos] - [results of neg])
nmfHmatAllByFact_combined_pos_negcoeff = nmfHmatAllByFact %>% 
  map(., ~as_tibble(.x, rownames = NA)) %>%
  map(., ~rownames_to_column(.x, "id")) %>%
  map(., ~pivot_longer(.x,
                       cols = contains("coeff"),
                       names_to = "factor",
                       values_to = "Weight")) %>%
  map(., ~mutate(.x, factor = gsub("_...coeff", "", factor))) %>% 
  map(., ~group_by(.x, factor, id)) %>%
  map(., ~summarise(.x, Weight = abs(.Primitive("-")(Weight[1], Weight[2])))) %>%
  map(., ~ungroup(.x)) %>%
  ## scale weights per signature so they add up to 1
  map(., ~group_by(.x, id)) %>% 
  map(., ~mutate(.x, sumWeight = sum(Weight))) %>% 
  map(., ~group_by(.x, id, factor)) %>% 
  map(., ~summarise(.x, Weight = Weight/sumWeight)) %>% 
  map(., ~ungroup(.x)) %>% 
  ## go to original format
  map(., ~pivot_wider(.x, names_from = factor, values_from = Weight)) %>%
  # arrange iterations (rows)
  map(., ~arrange(.x, id)) %>%
  map(., ~column_to_rownames(.x, 'id')) %>%
  # arrange dna repair factor names (columns)
  map(.f = list(. %>% select(coefficient_Resamp[[1]] %>% 
                               data.frame %>% 
                               select(contains("poscoeff")) %>% 
                               names() %>% 
                               ## when making it data.frame the '(' '>' ')' '×' symbols become '.', so revert this
                               #str_replace(., 
                               #           "([A,C,G,T])\\.([C,T])\\.([A,C,G,T])\\.([A,C,G,T])_", 
                               #           "\\1(\\2>\\3)\\4_") %>% 
                               #str_replace(., "\\.", "\\×") %>% 
                               gsub("_poscoeff", "", .)))) %>%
  map(., ~as.matrix(.x))



## Cluster the NMF factors, and find combination that yields best silhouette (and possibly deviation from random data)
## Also run in parallel on as many CPUs as specified
## Returns a data.table with nFact, k, and the quality measures per nFact*k combination.
#' @param maxNumClusters Global override - holds true for any number of nFact... note that some nFact will have less clusters than that 
#'                       For a given nFact, the max # clusters is 5*nFact.  Meaning: max 10 clusters for 2 factors, max 15 for 3 factors, etc.
#'                       (probably even this many does not make sense to extract...)

# initialize cosine distance function which is used for clustering
cosineDist <- function(x) {
  mat = as.dist(t(1 - cosine(t(x))))
  return(mat)
}

# initialize function for clustering
optimizeClustersPar = function(nmfHmatAllByFact_combined_pos_negcoeff = NULL, maxNumClusters = maxK, threads = nCPUs) {
  # the highest nFact (number of factors) we had checked in this dataset... the lowest is always assumed to be ==2 but this could be inspected as well
  maxNumFact = names(nmfHmatAllByFact_combined_pos_negcoeff) %>% 
    str_match("^nFact=(\\d+)") %>% .[, 2] %>% as.integer %>% max
  
  cl = makeCluster(threads, useXDR = F) # 'makeCluster' refers to computer clusters!
  clusterEvalQ(cl, library(cluster)) # to be able to run "pam" ('cluster' refers to data clustering, not computer clusters!)  
  sigDiscoScores = data.table()
  numNmfIters = -1 # this will be found from the nmfHmatAllByFact_combined_pos_negcoeff
  
  gc()
  
  for (aNmfSettings in names(nmfHmatAllByFact_combined_pos_negcoeff)) {
    
    # aNmfSettings is a string like "nFact=002", which encodes both the "nFact" parameter
    nFact = aNmfSettings %>% str_match("^nFact=(\\d+)") %>% .[, 2] %>% as.integer
    
    # how many iterations of NMF were run for this maxNumFact
    if (numNmfIters == -1) {
      numNmfIters = nmfHmatAllByFact_combined_pos_negcoeff[[aNmfSettings]] %>% rownames %>% str_match("^run(\\d+)_") %>% .[, 2] %>% as.integer %>% max
    } else {
      numNmfIters2 = nmfHmatAllByFact_combined_pos_negcoeff[[aNmfSettings]] %>% rownames %>% str_match("^run(\\d+)_") %>% .[, 2] %>% as.integer %>% max
      if (numNmfIters != numNmfIters2) {
        cat(sprintf("For aNmfSettings=%s, the numNmfIters %d is different than previously, when it was %d. Will continue running normally.\n", aNmfSettings, numNmfIters2, numNmfIters))
      }
    }
    
    ## get matrix of cosine distances between nFact×r factors (nFact = aNmfSettings, r = 1000 matrices)
    # each factor's vector goes from coord. origin to a point whose coordinates are its weights vector (so if 8 chromatin features, it's an intersection point of 8 dimensions)
    nmfReBootstrappedDists = cosineDist(nmfHmatAllByFact_combined_pos_negcoeff[[aNmfSettings]])
    gc()
    
    # parallel code: in each thread, runs the "pam" clustering for a different # clusters, starting from the same matrix of cosine similarities
    # (note: Sanger (Alexandrov et al.) always run this with #clusters = #NMF factors, which may not be ideal)
    clusterExport(cl = cl, list("nmfReBootstrappedDists"), envir = environment()) # , envir=tand)
    pamOutputByIter = parLapplyLB(# load balancing
      cl = cl,
      X = as.list(rev(2:min(maxNumClusters, nFact * 5))), # pam() with different values of k is run in parallel 
      fun = function(x) {
        set.seed(42 + x)
        # partitioning (clustering) of the nFact×r factors into k clusters "around medoids" (a more robust version of K-means), based on the cosine measures matrix
        pam(nmfReBootstrappedDists, k = x, diss = T)
      })
    gc()
    
    for (clusters in pamOutputByIter) {
      k = length(clusters$medoids)
      ## get clustering "stability": mean(clusters$silinfo$clus.avg.widths), i.e. average-across-K-clusters of the average silhouette widths per cluster
      # "silhouette width" of a sample (i.e. a given nFactor×run) i in a cluster is the "dissimilarity between i and the closest sample from an external cluster"
      sigDiscoScores = sigDiscoScores %>% 
        rbind(data.table(nFact = nFact, 
                         k = k,
                         avgClScore = mean(clusters$silinfo$clus.avg.widths), 
                         minClScore = min(clusters$silinfo$clus.avg.widths),
                         secondWorstClScore = clusters$silinfo$clus.avg.widths[order(clusters$silinfo$clus.avg.widths)][2],
                         avgPtScore = mean(silhouette(clusters)[, "sil_width"]), 
                         medianPtScore = median(silhouette(clusters)[, "sil_width"])))
      cat(sprintf("%02d\t%02d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                  nFact, k,
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
sigClusterScores = optimizeClustersPar(nmfHmatAllByFact_combined_pos_negcoeff, 
                                       maxNumClusters = maxK)

## Save clustering results 
saveRDS(sigClusterScores, paste0("heatmap_Ks/", "sigClusterScores.Rda")) # 'quality' of clustering(s), used to determine the number of signatures 
#sigClusterScores = readRDS(paste0("heatmap_Ks/", "sigClusterScores.Rda"))

## Visualize the clustering quality scores in heatmaps (facetted by the 'smooth' parameter)
heatmap_clustering = sigClusterScores[k <= maxK] %>% 
  melt(id.vars = c("nFact", "k"),
       measure.vars = "avgClScore") %>%
  rename("Stability" = "value")

heatmap_clustering_plot = ggplot(heatmap_clustering,
                                 aes(nFact, k)) +
  scale_x_continuous(breaks = seq(2, maxK)) +
  scale_y_continuous(breaks = seq(2, maxK)) +
  coord_fixed() +
  geom_tile(aes(fill = Stability)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0.4) +
  geom_label(aes(label = round(Stability, 2)),
             fontface = "bold",
             label.r = unit(0.3, "lines"),
             label.size = unit(0, "mm"),
             label.padding	= unit(0.15, "lines"),
             size = 5) +
  xlab("N factors") +
  ylab("K clusters") +
  theme_classic() +
  theme(text = element_text(size = 20))

ggsave("heatmap_Ks/NMF_heatmap_clustering.jpg",
       plot = heatmap_clustering_plot,
       device = "jpg",
       width = 12.5,
       height = 7,
       dpi = 900)


####
### Run NMF for optimal signatures

## for plots
dir.create("plots")

jet_colors = colorRampPalette(c("gray", "red", "yellow", "green", "cyan", "blue", "magenta", "black"))
## if you want a fixed assignment (i.e. that is constant across plots) of a given color to a given sample's feature:
levels = coefficient_table %>%
  ungroup %>% 
  separate(sample_id, into = c('type', 'sample_id'), sep = "_", fill = "warn") %>% 
  mutate(type = gsub("$", "\\-\\/\\-", type),
         type = gsub("normal.*", "wild-type", type)) %>% 
  select(type) %>% 
  distinct %>% 
  pull(type) 
  
fixed_jet_colors = jet_colors(length(levels))
names(fixed_jet_colors) = levels

# paint the weights with matching colors respect to the feature-/- colors in exposures
fixed_jet_colors_weights = fixed_jet_colors[-length(fixed_jet_colors)]
names(fixed_jet_colors_weights) = levels[-length(levels)] %>% gsub("\\-\\/\\-", "", .)

top_n_samples = 30


## create plots
                             # N,K
for(optimal_nFact_k in list(c(10,9))){
  
  ## extract K signatures (i.e, each cluster's medoid) for the desired combination of the number of NMF factors (nFact) (×1000 resampled matrices) and number of clusters (k)
  nFact = optimal_nFact_k[1]
  optimal_k = optimal_nFact_k[2] # the final number of signatures 

  wMatHere = nmfWmatAllByFact[[sprintf("nFact=%03d", nFact)]]
  hMatHere = nmfHmatAllByFact[[sprintf("nFact=%03d", nFact)]] # the NMF w-matrices had been transposed, so rows in the H-matrices and in the W-matrices are the same
  rownames(wMatHere) = gsub("run", "Run ", rownames(wMatHere))
  rownames(wMatHere) = gsub("fact0", "", rownames(wMatHere))
  rownames(wMatHere) = gsub("fact", "", rownames(wMatHere))
  rownames(wMatHere) = gsub("_", "\nFactor ", rownames(wMatHere))
  rownames(hMatHere) = gsub("run", "Run ", rownames(hMatHere))
  rownames(hMatHere) = gsub("fact0", "", rownames(hMatHere))
  rownames(hMatHere) = gsub("fact", "", rownames(hMatHere))
  rownames(hMatHere) = gsub("_", "\nFactor ", rownames(hMatHere))

  # need to run the clustering again...
  set.seed(42 + optimal_k) # this should reproduce the same output as in the "optimizeClustersPar"
  clustHere = pam(cosineDist(hMatHere), k = optimal_k, diss = T)
  # medoid: sample in a cluster that is the least dissimilar on average to all the other samples in that cluster
  signatures = copy(hMatHere[clustHere$medoids,])
  ## automatic ordering of signature names
  signature_names = sort(rownames(signatures))
  ## manual ordering of signature names
  signature_names = c("Run 049\nFactor 1",
                      "Run 079\nFactor 9",
                      "Run 066\nFactor 5",
                      "Run 047\nFactor 2",
                      "Run 094\nFactor 6",
                      "Run 027\nFactor 3",
                      "Run 065\nFactor 10",
                      "Run 082\nFactor 4",
                      "Run 055\nFactor 8")
  
  # get sample exposures for these cluster medoids (signatures)
  sample_exposures = wMatHere %>% 
    as_tibble() %>% 
    mutate(signature = rownames(wMatHere)) %>% 
    filter(signature %in% signature_names) %>% 
    arrange(signature) %>% relocate(signature)
  
  # and again, subtract results of pos - results of neg, and get the absolute...
  signatures_combined_pos_negcoeff = signatures %>% 
    as_tibble(rownames = NA) %>%
    rownames_to_column("id") %>%
    pivot_longer(cols = contains("coeff"),
                 names_to = "factor",
                 values_to = "Weight") %>%
    mutate(factor = gsub("_...coeff", "", factor)) %>% 
    group_by(factor, id) %>%
    summarise(Weight = abs(.Primitive("-")(Weight[1], Weight[2]))) %>%
    ungroup() %>%
    ## scale weights per signature so they add up to 1
    group_by(id) %>% 
    mutate(sumWeight = sum(Weight)) %>% 
    group_by(id, factor) %>% 
    summarise(Weight = Weight/sumWeight) %>% 
    ungroup() %>% 
    ## go to original format
    pivot_wider(names_from = factor, values_from = Weight) %>%
    # arrange iterations (rows)
    arrange(id) %>%
    column_to_rownames('id') %>%
    # arrange dna repair factor names (columns)
    select(coefficient_Resamp[[1]] %>% 
             data.frame %>% 
             select(contains("poscoeff")) %>% 
             names() %>% 
             gsub("_poscoeff", "", .)) %>%
    as.matrix()
  
  ## Signature exposures in samples
  # first write it for VAE
  dir.create("exposures")
  write_tsv(sample_exposures %>% 
              pivot_longer(cols = !contains("signature"), names_to = "Sample", values_to = "Exposure") %>% 
              rename("Signature" = "signature") %>% 
              filter(!is.na(Exposure)) %>% 
              mutate(Signature = gsub("\n| ", "_", Signature)),
            paste0("exposures/fct", nFact, "_", "k", optimal_k, "_exposures.tsv"))
  
  exposures = sample_exposures %>%
    pivot_longer(cols = !contains("signature"), names_to = "sample_id", values_to = "Exposure") %>% 
    mutate(sample_id = gsub("_nIter.*$", "", sample_id)) %>% 
    rename("Signature" = "signature") %>% 
    filter(!is.na(Exposure)) %>% 
    mutate(Signature = factor(Signature, levels = signature_names)) %>% 
    group_by(Signature) %>% 
    arrange(desc(Exposure)) %>% 
    slice_head(n = top_n_samples) %>% 
    ungroup %>% 
    mutate(sample_id = gsub("_.*", "\\-\\/\\-", sample_id),
           sample_id = gsub("normal.*", "wild-type", sample_id)) %>% 
    rename("Sample type" = "sample_id")
  
  # factor weights in signatures
  weights = signatures_combined_pos_negcoeff %>% 
    as_tibble() %>%
    mutate(Signature = rownames(signatures_combined_pos_negcoeff),
           Signature = factor(Signature, levels = signature_names)) %>% 
    pivot_longer(cols = !contains("Signature"), names_to = "factor", values_to = "weight") %>% 
    relocate(Signature) %>% 
    relocate(factor) %>% 
    rename("Feature" = "factor") %>% 
    mutate(`Feature` = factor(`Feature`,
                              levels = paste0("feature", seq(1:maxK)),
                              ordered = T)) %>% 
    arrange(`Feature`, Signature)
  

  ### plotting
  
  ## weights plot
  weights_plot = ggplot(weights %>%
                          # scale weights per signature so they add up to 1
                          group_by(Signature) %>% 
                          mutate(sumWeight = sum(weight)) %>% 
                          group_by(Signature, `Feature`) %>% 
                          summarise(Weight = weight/sumWeight) %>% 
                          ungroup, 
                        aes(x = Signature,
                            y = Weight)) +
    scale_y_continuous(expand = c(0, 0),
                       breaks = seq(0, 1, 0.25),
                       labels = function(x) sub("0+$", "", x)) +
    geom_col(aes(fill = `Feature`)) +
    scale_fill_manual(values = fixed_jet_colors_weights, drop = F) +
    guides(fill = guide_legend(override.aes = list(size=10))) +
    facet_grid(cols = vars(Signature), scales = "free") +
    theme_classic() +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          text = element_text(size = 25),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 20))
  
  ## exposures plot: for each nmf signature, plot only the top n samples with highest Exposure
  pos = position_jitter(w = 0.5, h = 0, seed = 1)
  exposures_plot = ggplot(exposures, 
                          aes(x = Signature,
                              y = Exposure)) +
    scale_y_log10() +
    geom_point(aes(fill = `Sample type`),
               size = 5,
               alpha = 0.7,
               shape = 21,
               position = pos) +
    scale_fill_manual(values = fixed_jet_colors, drop = F) +
    guides(fill = guide_legend(override.aes = list(size = 10, shape = 21))) +
    # geom_text_repel(aes(label = ),
    #                     size = 4,
    #                     force = 15,
    #                     position = pos,
    #                     max.overlaps = 1000000,
    #                     min.segment.length = 0.01) +
    facet_wrap(facets = vars(Signature), scales = "free", nrow = 1) +
    theme_classic() +
    xlab("") +
    ylab(bquote("% Exposure (" * log[10] * ", top-" * .(top_n_samples) ~ "samples)")) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_text(angle = 90, hjust = 0.5),
          text = element_text(size = 25),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.spacing = unit(12, "mm"),
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 20))
  
  combined_plots = plot_grid(NULL,
                             plot_grid(exposures_plot, NULL, nrow = 1, rel_widths = c(1, 0.04*(optimal_k/8))),
                             NULL,
                             plot_grid(weights_plot, NULL, nrow = 1, rel_widths = c(1, 0.04*(optimal_k/8))),
                             nrow = 4,
                             rel_heights = c(0.02, 1, -0.01, 0.7))
  ggsave(paste0("plots/NMF_exposures_weights_plot__nFact", nFact, "_K", optimal_k, ".jpeg"),
         plot = combined_plots,
         device = "jpeg",
         width = 21.3,
         height = 12,
         dpi = 900)
}
