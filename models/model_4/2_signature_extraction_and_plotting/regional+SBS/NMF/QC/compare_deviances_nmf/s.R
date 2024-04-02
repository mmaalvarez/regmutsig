library(tidyverse)
library(data.table)
library(parallel)
library(msm) # rtnorm
#install.packages('terra', repos='https://rspatial.r-universe.dev')
library(NMF) # for NMF to determine optimal k signatures
#ENMTools::install.extras("NMF")
library(lsa)
#library(cluster)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("cosine", "lsa")
conflict_prefer("clusterExport", "parallel")
conflict_prefer("clusterEvalQ", "parallel")
conflict_prefer("parLapply", "parallel")
conflict_prefer("parLapplyLB", "parallel")
conflict_prefer("nmf", "NMF")


##### samples info

K562 = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/marcel_K562/metadata/WGS_clones_info.tsv") %>% 
  rename("Sample" = "sample_id",
         "alteration" = "genotype KOs") %>%
  mutate(dataset = "Marcel",
         Sample = gsub("_REP1", "", Sample),
         # I consider MSH3-/- as MMRwt
         altered_pathway_or_treatment_type = ifelse((str_detect(alteration, "MSH2|MSH6|MLH1|PMS2") | str_detect(`MMR deficiency expected`, "strong|middle")) & str_detect(alteration, "OGG1|MUTYH|NTH1|NEIL2"),
                                                    "MMR & BER",
                                                    ifelse(str_detect(alteration, "MSH2|MSH6|MLH1|PMS2") | str_detect(`MMR deficiency expected`, "strong|middle"),
                                                           "MMR",
                                                           ifelse(str_detect(alteration, "OGG1|MUTYH|NTH1|NEIL2"),
                                                                  "BER",
                                                                  "control"))),
         alteration = ifelse(alteration != "WT",
                             gsub("$", " KO", alteration),
                             alteration)) %>%  
  select(Sample, dataset, alteration, altered_pathway_or_treatment_type)

iPSC = c("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/kucab_2019/processed/sample_treatments.tsv",
         "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/zou_2021/processed/sample_gene_ko.tsv") %>% 
  # only Sample and info* columns are selected
  map_df(~read_tsv(.x)) %>% 
  rename("alteration" = "info1",
         "altered_pathway_or_treatment_type" = "info2",
         "Sample" = "sample_id") %>% 
  mutate(altered_pathway_or_treatment_type = gsub("^[a-k]_", "", altered_pathway_or_treatment_type),
         altered_pathway_or_treatment_type = gsub("Control", "control", altered_pathway_or_treatment_type),
         dataset = ifelse(str_detect(Sample, "MSM0"),
                          "Kucab",
                          ifelse(str_detect(Sample, "MSK0"),
                                 "Zou",
                                 "ERROR: Unexpected sample name")),
         altered_pathway_or_treatment_type = gsub("DNA damage response inhibitors", "DNA damage resp. inh.", altered_pathway_or_treatment_type),
         # "MMRko" sample with very low mut burden, treat as control
         `altered_pathway_or_treatment_type` = ifelse(Sample == "MSK0.123_s1",
                                                      "control",
                                                      `altered_pathway_or_treatment_type`),
         alteration = ifelse(dataset == "Zou",
                             gsub("$", " KO", alteration),
                             alteration)) %>% 
  select(Sample, dataset, alteration, altered_pathway_or_treatment_type)

petljak = read_tsv("/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/cell_lines/petljak_2022/info/metadata.tsv") %>% 
  rename("Sample" = "sample_id",
         "alteration" = "info1",
         "cell_line" = "info2") %>%
  mutate(dataset = "Petljak",
         altered_pathway_or_treatment_type = ifelse(alteration == "WT",
                                                    "control",
                                                    "APOBEC"),
         alteration = gsub("_KO", " KO", alteration)) %>% 
  select(Sample, dataset, alteration, altered_pathway_or_treatment_type)

# merge datasets metadata
samples_info = bind_rows(K562, iPSC) %>% 
  bind_rows(petljak)



## load res
critical_value = 1.28 # for CI80%
# 1.96 for CI95%

coefficient_table = lapply(c("../../../../1_parser_and_regressions/res/results_regional_feature_regressions_all_samples.tsv",
                             "../../../../../model_2/1_parser_and_regressions/res/results_SBS96_regressions_all_samples.tsv"),
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
  left_join(samples_info) %>% 
  select(Sample, dataset, alteration, altered_pathway_or_treatment_type, feature_name, non_ref_bin_feature_level, feature_type, stat, val) %>%
  pivot_wider(names_from = stat, values_from = val) %>% 
  arrange(feature_name, Sample) %>% 
  group_by(feature_name, Sample)
gc()


#### permutations

## Parameters and initializing of some objects
totalNumIters = 100
coefficient_Resamp = list()
set.seed(1)

## Generate matrices resampling betas from their CI % distributions (instead of UPmultinomial)
resample_from_CI = function(coefficient_table){
  
  ## for QC, the regional features will be completely random values (i.e. meaningless variables)
  coefficient_table = ungroup(coefficient_table) %>% 
    group_by(feature_type)
  
  coefficient_table_split = coefficient_table %>% 
    group_split()
  
  names(coefficient_table_split) = group_keys(coefficient_table)$feature_type
  
  # do actual sampling from CI in SBS96...
  coefficient_table_SBS96 = coefficient_table_split$SBS96 %>% 
    group_by(feature_name, Sample) %>% 
    summarise(resampled_estimate = rtnorm(n = 1,
                                          mean = estimate,
                                          sd = 1, 
                                          lower = conf.low,
                                          upper = conf.high))
  # ...but shuffle Sample__feature_name, so `estimate conf.low conf.high` row-wise for regional features are shuffled
  coefficient_table_regfeat = coefficient_table_split$regional_feature %>% 
    select(Sample, feature_name, estimate, contains("conf.")) %>% 
    unite("id", Sample, feature_name, sep = "__") %>% 
    transform(id = sample(id)) %>% 
    separate(id, into = c("Sample", "feature_name"), sep = "__") %>% 
    relocate(feature_name) %>% relocate(Sample) %>% 
    group_by(feature_name, Sample) %>% 
    summarise(resampled_estimate = rtnorm(n = 1,
                                          mean = estimate,
                                          sd = 1, 
                                          lower = conf.low,
                                          upper = conf.high))
  
  # return the 2 subtables rebound together
  bind_rows(coefficient_table_SBS96,
            coefficient_table_regfeat)
}

for (nIter in 1:totalNumIters) {
  print(paste0("nIter ", nIter))
  # for each sample (row) resample coefficients from CI % distrs.
  coefficient_Resamp[[nIter]] = resample_from_CI(coefficient_table) %>% 
    mutate("nIter" = nIter)
  gc()
}

# bind all permuted tables
coefficient_Resamp = bind_rows(coefficient_Resamp)

dir.create("heatmap_Ks")
write_tsv(coefficient_Resamp, paste0("heatmap_Ks/permuted_coefficients_", totalNumIters, "iters.tsv"))
#coefficient_Resamp = read_tsv("heatmap_Ks/permuted_coefficients_100iters.tsv")



#####################################################################
#### NMF
#####################################################################


#### evaluate best combination of n variables and k signatures

coefficient_Resamp = coefficient_Resamp %>% 
  ungroup %>% 
  pivot_wider(names_from = feature_name, values_from = resampled_estimate) %>% 
  unite("Sample", Sample, nIter, sep = "__") %>% 
  # NAs to 0
  mutate_all(~ifelse(is.na(.), 0, .)) %>% 
  # remove features with all 0s
  select(where(~any(. != 0)))

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
                           by = "Sample",
                           suffixes = c("_poscoeff", "_negcoeff")) %>%
  mutate(Sample = gsub("__muts", "_muts", Sample)) %>% 
  separate(Sample, into = c("id", "nIter"), sep = "__", remove = F) %>% 
  mutate(Sample = gsub("__", "_nIter", Sample)) %>% 
  column_to_rownames("Sample") %>% 
  select(-id) %>% 
  split(., f = .$nIter) %>% 
  map(.f = list(. %>% data.matrix))



## Run NMF for each resampled coefficients matrix generated in the previous step

## Parameters and initializing of some objects
set.seed(1)
maxK = 12 # max number of signatures to consider, it will go from 2 to maxK -- shouldn't be larger than nº of features + expected_n_cosmic
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
  
  # ...from each of the 1000 permuted matrices
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
                               ### NMF
                               nmf_run = NMF::nmf(x_nozerocols, 
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
                               # when making it data.frame the '(' '>' ')' '×' symbols become '.', so revert this
                               full_matrix_colnames = colnames(full_matrix) %>% 
                                 str_replace(., 
                                             "([A,C,G,T])\\.([C,T])\\.([A,C,G,T])\\.([A,C,G,T])_", 
                                             "\\1(\\2>\\3)\\4_") %>% 
                                 str_replace(., "X53BP1", "53BP1")
                               colnames(full_matrix) = full_matrix_colnames
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
#rm(idString, nmfOutput, nFact)

## Save clustering results
saveRDS(nmfHmatAllByFact, paste0("heatmap_Ks/", "Hmatrix.Rda")) # signatures
#nmfHmatAllByFact = readRDS(paste0("heatmap_Ks/", "Hmatrix.Rda"))
saveRDS(nmfWmatAllByFact, paste0("heatmap_Ks/", "Wmatrix.Rda")) # exposures
#nmfWmatAllByFact = readRDS(paste0("heatmap_Ks/", "Wmatrix.Rda"))


## [results of pos] + -[results of neg] to have a "-" sign when the larger/non-zero value comes from the neg. submatrix
nmfHmatAllByFact_combined_pos_negcoeff = nmfHmatAllByFact %>% 
  map(., ~as_tibble(.x, rownames = NA)) %>%
  map(., ~rownames_to_column(.x, "id")) %>%
  map(., ~pivot_longer(.x,
                       cols = contains("coeff"),
                       names_to = "dna_repair_feature_name",
                       values_to = "Weight")) %>%
  map(., ~mutate(.x, Weight = ifelse(str_detect(dna_repair_feature_name, "_negcoeff"),
                                     -Weight,
                                     Weight))) %>% 
  map(., ~mutate(.x, dna_repair_feature_name = gsub("_...coeff", "", dna_repair_feature_name))) %>% 
  map(., ~group_by(.x, dna_repair_feature_name, id)) %>%
  map(., ~summarise(.x, Weight = .Primitive("+")(Weight[1], Weight[2]))) %>%
  map(., ~ungroup(.x)) %>%
  ## scale weights per signature so they add up to 1
  #map(., ~group_by(.x, id)) %>% 
  #map(., ~mutate(.x, sumWeight = sum(Weight))) %>% 
  #map(., ~group_by(.x, id, dna_repair_feature_name)) %>% 
  #map(., ~summarise(.x, Weight = Weight/sumWeight)) %>% 
  #map(., ~ungroup(.x)) %>% 
  ## go to original format
  map(., ~pivot_wider(.x, names_from = dna_repair_feature_name, values_from = Weight)) %>%
  # arrange iterations (rows)
  map(., ~arrange(.x, id)) %>%
  map(., ~column_to_rownames(.x, 'id')) %>%
  # arrange dna repair feature_name names (columns)
  map(.f = list(. %>% select(coefficient_Resamp[[1]] %>% 
                               data.frame %>% 
                               select(contains("poscoeff")) %>% 
                               names() %>% 
                               # when making it data.frame the '(' '>' ')'symbols become '.', so revert this
                               str_replace(., 
                                           "([A,C,G,T])\\.([C,T])\\.([A,C,G,T])\\.([A,C,G,T])_", 
                                           "\\1(\\2>\\3)\\4_") %>% 
                               str_replace(., "X53BP1", "53BP1") %>% 
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

# cosine SIMILARITY will be used later for similarities between signatures' regional features weight profiles
cosineSimil <- function(x) {
  mat = t(cosine(t(x)))
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
saveRDS(sigClusterScores, paste0("heatmap_Ks/", "sigClusterScores.Rda")) # 'quality' of clustering(s), used to defeature_nameine the number of signatures 
#sigClusterScores = readRDS(paste0("heatmap_Ks/", "sigClusterScores.Rda"))


## Visualize the clustering quality scores in heatmaps (faceted by the 'smooth' parameter)
heatmap_clustering = sigClusterScores[k <= maxK] %>% 
  melt(id.vars = c("nFact", "k"),
       measure.vars = "avgClScore") %>%
  rename("Stability" = "value")

heatmap_clustering_plot = ggplot(heatmap_clustering %>% 
                                   # only 10 factors and k
                                   filter(nFact<=maxK & k<=maxK),
                                 aes(nFact, k)) +
  scale_x_continuous(breaks = seq(2, maxK)) +
  scale_y_continuous(breaks = seq(2, maxK)) +
  #coord_fixed() +
  geom_tile(aes(fill = Stability)) +
  scale_fill_gradient2(low = "white", high = "blue", midpoint = 0,
                       breaks = seq(-0.25,1.25,0.25)) +
  geom_label(aes(label = round(Stability, 2)),
             #fontface = "bold",
             label.r = unit(0.01, "lines"),
             label.size = unit(0, "mm"),
             label.padding  = unit(0.01, "lines"),
             size = 8) +
  xlab("N factors") +
  ylab("K clusters") +
  theme_classic() +
  theme(text = element_text(size = 28),
        legend.position = "top",
        legend.text = element_text(angle = 90, vjust = 0.5, hjust = 0.25))

ggsave("heatmap_Ks/NMF_heatmap_clustering_maxnFact10_maxK10.jpg",
       plot = heatmap_clustering_plot,
       device = "jpg",
       width = 12.5,
       height = 7,
       dpi = 600)


####
### Run NMF for optimal signatures

## for QC (the deviance between NMF fitted model and original coefficients matrices in regfeat+SBS96 matrices vs. in SBS96-only matrices) AND shuffled reg feat coefficients (here)
deviances_list = list()
# this stores deviances after removing the reg feat reasiduals, for a direct comparison to SNS96only
deviances_NO_REG_FEAT_list = list()

# define combinations of nFact and k
range_nFact = seq(3, maxK)

for(optimal_nFact in range_nFact){
  
  nFact = optimal_nFact

  wMatHere = nmfWmatAllByFact[[sprintf("nFact=%03d", nFact)]]
  hMatHere = nmfHmatAllByFact[[sprintf("nFact=%03d", nFact)]] # the NMF w-matrices had been transposed, so rows in the H-matrices and in the W-matrices are the same
  
  #### for QC: calculate the deviance between each NMF fitted model and its corresponding original (positivized) coefficient matrix (i.e. for each nFact==[2-12] × iter==[1-100] combination); this will be compared to the equivalent when there are only SBS96 coeff, to see whether we actually get more accurate NMF fit when SBS96 coefficients are combined with regional feature coefficients
  deviances_list[[paste0('nFact',nFact)]] = list()
  deviances_NO_REG_FEAT_list[[paste0('nFact',nFact)]] = list()
  
  for(nIter in 1:totalNumIters){
    
    original_matrix = coefficient_Resamp[[nIter]]
    # remove first column
    original_matrix = original_matrix[, colnames(original_matrix) != "nIter"]
    # sometimes it can be that e.g. coefficient_Resamp[[6]] contains ..._nIter013 in rownames
    actual_nIter = rownames(original_matrix) %>% gsub(".*_nIter", "", .) %>% unique %>% as.numeric %>% sprintf("%03d", .)
    
    exposures_this_actual_nIter = wMatHere %>% 
      data.frame %>% 
      rownames_to_column("nIter_nFact") %>% 
      filter(str_detect(nIter_nFact, paste0("run",actual_nIter,"_"))) %>% 
      column_to_rownames("nIter_nFact") %>% 
      as.matrix %>% 
      t
    
    weights_this_actual_nIter = hMatHere %>% 
      data.frame %>% 
      rownames_to_column("nIter_nFact") %>% 
      filter(str_detect(nIter_nFact, paste0("run",actual_nIter,"_"))) %>% 
      column_to_rownames("nIter_nFact") %>% 
      as.matrix
    
    exposures_by_weights = exposures_this_actual_nIter %*% weights_this_actual_nIter
    
    # mean of the squared differences
    original_vs_nmf_deviance = sum(as.numeric((exposures_by_weights - original_matrix)^2)) / prod(dim(original_matrix))
    
    deviances_list[[paste0('nFact',nFact)]][[paste0('nIter',nIter)]] = original_vs_nmf_deviance
    
    # mean of the squared differences without regfeat columns
    original_matrix_NO_REG_FEAT = data.frame(original_matrix) %>% 
      select(contains(".")) %>% as.matrix
    exposures_by_weights_NO_REG_FEAT = data.frame(exposures_by_weights) %>% 
      select(contains(".")) %>% as.matrix
    original_vs_nmf_deviance_NO_REG_FEAT = sum(as.numeric((exposures_by_weights_NO_REG_FEAT - original_matrix_NO_REG_FEAT)^2)) / prod(dim(original_matrix_NO_REG_FEAT))
    deviances_NO_REG_FEAT_list[[paste0('nFact',nFact)]][[paste0('nIter',nIter)]] = original_vs_nmf_deviance_NO_REG_FEAT
    
    gc()
  }
  gc()
}

deviances_with_shuffled_regfeat_coeff = bind_rows(deviances_list) %>% 
  as.matrix %>% t %>% data.frame %>% `colnames<-`(paste0("nFact", range_nFact))
#deviances_with_shuffled_regfeat_coeff = read_tsv("deviances_with_shuffled_regfeat_coeff.tsv")
write_tsv(deviances_with_shuffled_regfeat_coeff, 
          "deviances_with_shuffled_regfeat_coeff.tsv")
deviances_with_shuffled_regfeat_coeff = mutate(deviances_with_shuffled_regfeat_coeff,
                                               type = "Regional features (shuffled) + SBS96")

deviances_with_NO_shuffled_REG_FEAT = bind_rows(deviances_NO_REG_FEAT_list) %>% 
  as.matrix %>% t %>% data.frame %>% `colnames<-`(paste0("nFact", range_nFact))
#deviances_with_NO_shuffled_REG_FEAT = read_tsv("deviances_with_shuffled_REG_FEAT_REMOVED.tsv")
write_tsv(deviances_with_NO_shuffled_REG_FEAT, 
          "deviances_with_shuffled_REG_FEAT_REMOVED.tsv")
deviances_with_NO_shuffled_REG_FEAT = mutate(deviances_with_NO_shuffled_REG_FEAT,
                                             type = "Regional features (shuffled) REMOVED + SBS96")

# load SBS96-only NMF deviances
deviances_without_regfeat_coeff = read_tsv("../../../SBS/NMF/QC/deviances_without_regfeat_coeff.tsv") %>% 
  mutate(type = "SBS96 only")

# load Reg feat (no shuffle) + SBS96 NMF deviances
deviances_with_regfeat_coeff = read_tsv("deviances_with_regfeat_coeff.tsv") %>% 
  mutate(type = "Regional features + SBS96")

# load Reg feat (no shuffle) REMOVED + SBS96 NMF deviances
deviances_with_regfeat_coeff_REMOVED = read_tsv("deviances_REGFEAT_REMOVED.tsv") %>% 
  mutate(type = "Regional features REMOVED + SBS96")

### comparison of deviances
deviances_with_without_regfeat_coeff = deviances_with_regfeat_coeff %>% # Reg feat + SBS96
  bind_rows(deviances_with_regfeat_coeff_REMOVED,                       # (Reg feat REMOVED) SBS96
            deviances_with_shuffled_regfeat_coeff,                      # Reg feat (shuffled) + SBS96
            deviances_with_NO_shuffled_REG_FEAT,                        # (Reg feat shuffled REMOVED) SBS96
            deviances_without_regfeat_coeff) %>%                        # SBS96 only
  mutate(type = factor(type, levels = c("Regional features + SBS96",
                                        "Regional features (shuffled) + SBS96",
                                        "Regional features REMOVED + SBS96",
                                        "Regional features (shuffled) REMOVED + SBS96",
                                        "SBS96 only"), ordered = T))

plot_comparison_deviances = deviances_with_without_regfeat_coeff %>% 
  pivot_longer(cols = !matches("type"), names_to = "nFact", values_to = "deviance") %>% 
  mutate(nFact = factor(nFact, levels = paste0("nFact", range_nFact), ordered = T)) %>% 
  ggplot(aes(x = type,
             y = deviance)) +
  geom_violin(aes(fill = type)) +
  geom_boxplot(width = 0.3) +
  scale_fill_manual(values = c("red", "lightblue", "orange", "turquoise", "yellow")) +
  ## mann whitney between all samples
  ggsignif::geom_signif(comparisons = c(combn(c("Regional features + SBS96",
                                              "Regional features (shuffled) + SBS96"), 2, simplify = F),
                                        combn(c("Regional features REMOVED + SBS96",
                                                "Regional features (shuffled) REMOVED + SBS96",
                                                "SBS96 only"), 2, simplify = F)),
                        step_increase = 0.05,
                        map_signif_level = T, # asterisks if TRUE
                        test = "wilcox.test",
                        test.args = "two.sided") +  # or "less", "greater"
  facet_wrap(facets = vars(nFact),
             nrow = 1,
             strip.position="bottom") +
  theme_bw() +
  xlab("") +
  ylab("Σ (i.coefficients − i.exposures × i.weights)² / (nSamples × mFeatures) ; i ∈ [1,100]") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "top")
ggsave(plot = plot_comparison_deviances, 
       "comparison_deviances.jpg",
       device = "jpeg",
       width = 14.2,
       height = 8,
       dpi = 600)
