library(tidyverse)
library(data.table)
library(parallel)
library(msm) # rtnorm
#install.packages('terra', repos='https://rspatial.r-universe.dev')
library(NMF) # for NMF to defeature_nameine optimal k signatures
#ENMTools::install.extras("NMF")
library(lsa)
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
         # I consider MSH3-/- as MMRwt
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

nmfHmatAllByFact_halves = list()
nmfWmatAllByFact_halves = list()

args = commandArgs(trailingOnly=TRUE)
genome_half = ifelse(interactive(),
                     yes = "1st_half", #"2nd_half", "shifted")
                     no = args[1])

coefficient_table = lapply(c(paste0("/g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model_4/2_signature_extraction_and_plotting/regional+SBS/NMF/QC/consistency/", genome_half, "_genome/res/results_regional_feature_regressions_all_samples.tsv"),
                             paste0("/g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/models/model_4/2_signature_extraction_and_plotting/regional+SBS/NMF/QC/consistency/", genome_half, "_genome/res/results_SBS96_regressions_all_samples.tsv")),
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
  group_by(Sample, feature_name)
gc()


#### permutations

## Parameters and initializing of some objects
totalNumIters = 100
coefficient_Resamp = list()
set.seed(1)

## Generate matrices resampling betas from their CI % distributions (instead of UPmultinomial)
resample_from_CI = function(coefficient_table){
  summarise(coefficient_table, resampled_estimate = rtnorm(n = 1,
                                                           mean = estimate,
                                                           sd = 1, 
                                                           lower = conf.low,
                                                           upper = conf.high))}

if(identical(Sys.glob(paste0("heatmap_Ks/permuted_coefficients_", totalNumIters, "iters_", genome_half,".tsv")), character(0))){
  
  ## no permuted coefficients are done yet, so do them
  
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
  write_tsv(coefficient_Resamp, paste0("heatmap_Ks/permuted_coefficients_", totalNumIters, "iters_", genome_half,".tsv"))
  
} else {
  # load it
  coefficient_Resamp = read_tsv(paste0("heatmap_Ks/permuted_coefficients_", totalNumIters, "iters_", genome_half,".tsv"))
}


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
    # WARNING: the rbind causes that the colnames become always ..._nIter1, e.g. MSM0.96_nIter2 --> MSM0.96_nIter1; this shouldn't matter as long as all samples are always in the same column i
    nmfWmatAllByFact[[idString]] = nmfWmatAllByFact[[idString]] %>% rbind(t(nmfOutput@fit@W)) # same thing that you get using the basis(NMFfit) command in R  -> transformed data
    gc()
  }
  gc()
}
#stopCluster(cl)
#rm(idString, nmfOutput, nFact)

## Save clustering results
saveRDS(nmfHmatAllByFact, paste0("heatmap_Ks/", "Hmatrix_", genome_half, ".Rda")) # signatures
#nmfHmatAllByFact = readRDS(paste0("heatmap_Ks/", "Hmatrix_", genome_half, ".Rda"))
nmfHmatAllByFact_halves[[genome_half]] = nmfHmatAllByFact

saveRDS(nmfWmatAllByFact, paste0("heatmap_Ks/", "Wmatrix_", genome_half, ".Rda")) # exposures
#nmfWmatAllByFact = readRDS(paste0("heatmap_Ks/", "Wmatrix_", genome_half, ".Rda"))
nmfWmatAllByFact_halves[[genome_half]] = nmfWmatAllByFact


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
saveRDS(sigClusterScores, paste0("heatmap_Ks/", "sigClusterScores_", genome_half, ".Rda")) # 'quality' of clustering(s), used to defeature_nameine the number of signatures 
#sigClusterScores = readRDS(paste0("heatmap_Ks/", "sigClusterScores_", genome_half, ".Rda"))


## Visualize the clustering quality scores in heatmaps (faceted by the 'smooth' parameter)
heatmap_clustering = sigClusterScores[k <= maxK] %>% 
  melt(id.vars = c("nFact", "k"),
       measure.vars = "avgClScore") %>%
  rename("Stability" = "value")

heatmap_clustering_plot = ggplot(heatmap_clustering %>% 
                                   # only maxK factors and k
                                   filter(nFact<=maxK & k<=maxK),
                                 aes(nFact, k)) +
  scale_x_continuous(breaspearman = seq(2, maxK)) +
  scale_y_continuous(breaspearman = seq(2, maxK)) +
  #coord_fixed() +
  geom_tile(aes(fill = Stability)) +
  scale_fill_gradient2(low = "white", high = "blue", midpoint = 0,
                       breaspearman = seq(-0.25,1.25,0.25)) +
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

ggsave(paste0("heatmap_Ks/NMF_heatmap_clustering_maxnFact", maxK, "_maxK", maxK, "_", genome_half, ".jpg"),
       plot = heatmap_clustering_plot,
       device = "jpg",
       width = 12.5,
       height = 7,
       dpi = 600)



#### compare NMF results between both genome halves AND null data


## load all results

#coefficient_Resamp_halves = list()
#nmfHmatAllByFact_halves = list()
#nmfWmatAllByFact_halves = list()

for(genome_half in c("1st_half", "2nd_half", "null")){
  
  coefficient_Resamp_halves[[genome_half]] = read_tsv(paste0("heatmap_Ks/permuted_coefficients_100iters_", genome_half, ".tsv"))
  nmfHmatAllByFact_halves[[genome_half]] = readRDS(paste0("heatmap_Ks/Hmatrix_", genome_half, ".Rda"))
  nmfWmatAllByFact_halves[[genome_half]] = readRDS(paste0("heatmap_Ks/Wmatrix_", genome_half, ".Rda"))
}


### compare sample exposures of a)1st vs 2nd, and b) 1st vs null, by nFact and run
nmfWmatAllByFact_halves_bound_list = list()
for(genome_half in names(nmfWmatAllByFact_halves)){
  nmfWmatAllByFact_halves_bound_list[[genome_half]] = nmfWmatAllByFact_halves[[genome_half]] %>% 
    map2(names(nmfWmatAllByFact_halves[[genome_half]]),
         ~mutate(data.frame(.x), nFact = .y)) %>% 
    bind_rows %>% 
    relocate(nFact) %>% 
    rownames_to_column("run_fact") %>% 
    mutate(run_fact = gsub("\\.\\.\\..*", "", run_fact)) %>% 
    separate(run_fact, into = c("run", "fact"), sep = "_")
}
nmfWmatAllByFact_halves_bound = nmfWmatAllByFact_halves_bound_list %>% 
  map2(names(nmfWmatAllByFact_halves),
       ~mutate(data.frame(.x), genome_half = .y)) %>% 
  bind_rows %>% 
  relocate(genome_half) %>% 
  pivot_longer(cols = -c(genome_half, run, fact, nFact), names_to = "Sample", values_to = "Exposure") %>% 
  pivot_wider(names_from = genome_half, values_from = Exposure) %>%
  mutate(nFact = factor(nFact, levels = paste0("nFact=", str_pad(seq(2,maxK), 3, pad = "0"))),
         fact = factor(fact, levels = paste0("fact", str_pad(seq(maxK,1), 2, pad = "0"))))

### WARNING: in 2 NMF (diff. seed or slightly diff data) for a given nFact, the factor called "factor01" always refers to the same factor, etc.
# test_compare_same_factor_diff_runs = nmfWmatAllByFact_halves_bound_list %>% 
#   map2(names(nmfWmatAllByFact_halves),
#        ~mutate(data.frame(.x), genome_half = .y)) %>% 
#   bind_rows %>% 
#   relocate(genome_half) %>% 
#   pivot_longer(cols = -c(genome_half, run, fact, nFact), names_to = "Sample", values_to = "Exposure") %>% 
#   pivot_wider(names_from = run, values_from = Exposure) %>%
#   mutate(nFact = factor(nFact, levels = paste0("nFact=", str_pad(seq(2,maxK), 3, pad = "0"))),
#          fact = factor(fact, levels = paste0("fact", str_pad(seq(maxK,1), 2, pad = "0")))) %>% 
#   filter(genome_half == "1st_half") %>% 
#   select(nFact, fact, run001, run002)
# ggplot(test_compare_same_factor_diff_runs,
#        aes(x = `run001`,
#            y = `run002`)) +
#   geom_point(alpha = 0.5,
#              size = 0.5) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   facet_grid(rows = vars(fact),
#              cols = vars(nFact),
#              scales = "free",
#              space = "free",
#              switch = "both")


### DO THIS for a single run within each genome_half/nFact, as in diff. runs a factor can be called fact01 in one run and fact02 in another run, e.g.

### AND slice_max the exposure in BOTH first_half and 2nd_half, to make sure of this

max_exp_fact_1st_half = nmfWmatAllByFact_halves_bound %>%
  ## run001 only
  filter(run == "run001") %>% 
  select(nFact, Sample, `1st_half`) %>% 
  group_by(nFact, Sample) %>% 
  slice_max(`1st_half`) %>% 
  ungroup
  
max_exp_fact_2nd_half = nmfWmatAllByFact_halves_bound %>%
  ## run001 only
  filter(run == "run001") %>% 
  select(nFact, Sample, `2nd_half`) %>% 
  group_by(nFact, Sample) %>% 
  slice_max(`2nd_half`) %>% 
  ungroup

max_exp_fact_null = nmfWmatAllByFact_halves_bound %>%
  ## run001 only
  filter(run == "run001") %>% 
  select(nFact, Sample, `null`) %>% 
  group_by(nFact, Sample) %>% 
  slice_max(`null`) %>% 
  ungroup

max_exp_fact_all = left_join(max_exp_fact_1st_half,
                                  max_exp_fact_2nd_half) %>% 
  left_join(max_exp_fact_null)

# spearman corr
exposures_plot_table_spearman_1st_vs_2nd = max_exp_fact_all %>% 
  split(as.factor(.$nFact)) %>%
  map(~cor(.x$`1st_half`, .x$`2nd_half`, method = "spearman")) %>%
  map(~broom::tidy(.x)) %>%
  map(~pull(.x, x)) %>%
  unlist %>%
  data.frame %>%
  rownames_to_column("nFact") %>%
  rename("Odd vs. even chr\nSpearman corr" = ".") %>%
  mutate(nFact = factor(nFact, levels = paste0("nFact=", str_pad(seq(2,maxK), 3, pad = "0"))))

exposures_plot_table_spearman_1st_vs_null = max_exp_fact_all %>%
  split(as.factor(.$nFact)) %>%
  map(~cor(.x$`1st_half`, .x$`null`, method = "spearman")) %>%
  map(~broom::tidy(.x)) %>%
  map(~pull(.x, x)) %>%
  unlist %>%
  data.frame %>%
  rownames_to_column("nFact") %>%
  rename("Odd chr vs. null\nSpearman corr" = ".") %>%
  mutate(nFact = factor(nFact, levels = paste0("nFact=", str_pad(seq(2,maxK), 3, pad = "0"))))

# plot odd vs even chr
plot_exposures_odd_vs_even_chr = ggplot(left_join(max_exp_fact_all,
                                                  exposures_plot_table_spearman_1st_vs_2nd),
                                        aes(x = `1st_half`,
                                            y = `2nd_half`,
                                            col = `Odd vs. even chr\nSpearman corr`)) +
  geom_point(alpha = 0.5,
             size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_tile() +
  scale_color_gradient(low = "lightblue",
                       high = "blue") +
  geom_label(data = exposures_plot_table_spearman_1st_vs_2nd,
             aes(label = round(`Odd vs. even chr\nSpearman corr`, 2),
                 x = 0.8, y = 0.2),
             col = "black",
             label.padding = unit(0.1, "lines"),
             label.size = unit(0, "lines"),
             size = 2) +
  facet_wrap(facets = vars(nFact),
             nrow = 1) +
  xlab("NMF (single run) exposures of the highest-exposure factor, using the ODD chromosomes") +
  ylab("NMF (single run) exposures of the highest-exposure factor, using the EVEN chromosomes") +
  theme_classic() +
  theme(text = element_text(size = 6))
ggsave("plot_exposures_odd_vs_even_chr.jpg",
       plot = plot_exposures_odd_vs_even_chr,
       device = "jpg",
       width = 16,
       height = 9,
       dpi = 600)

# odd vs. null
plot_exposures_odd_vs_null = ggplot(left_join(max_exp_fact_all,
                                              exposures_plot_table_spearman_1st_vs_null),
                                        aes(x = `1st_half`,
                                            y = `null`,
                                            col = `Odd chr vs. null\nSpearman corr`)) +
  geom_point(alpha = 0.5,
             size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_tile() +
  scale_color_gradient(low = "lightblue",
                       high = "blue") +
  geom_label(data = exposures_plot_table_spearman_1st_vs_null,
             aes(label = round(`Odd chr vs. null\nSpearman corr`, 2),
                 x = 0.8, y = 0.2),
             col = "black",
             label.padding = unit(0.1, "lines"),
             label.size = unit(0, "lines"),
             size = 2) +
  facet_wrap(facets = vars(nFact),
             nrow = 1) +
  xlab("NMF (single run) exposures of the highest-exposure factor, using the ODD chromosomes") +
  ylab("NMF (single run) exposures of the highest-exposure factor, using the NULL data") +
  theme_classic() +
  theme(text = element_text(size = 6))
ggsave("plot_exposures_odd_vs_null.jpg",
       plot = plot_exposures_odd_vs_null,
       device = "jpg",
       width = 16,
       height = 9,
       dpi = 600)


###############################################################################

# trinuc_96 = c("A(C>A)A", "A(C>A)C", "A(C>A)G", "A(C>A)T",
#               "C(C>A)A", "C(C>A)C", "C(C>A)G", "C(C>A)T",
#               "G(C>A)A", "G(C>A)C", "G(C>A)G", "G(C>A)T",
#               "T(C>A)A", "T(C>A)C", "T(C>A)G", "T(C>A)T",
#               "A(C>G)A", "A(C>G)C", "A(C>G)G", "A(C>G)T",
#               "C(C>G)A", "C(C>G)C", "C(C>G)G", "C(C>G)T",
#               "G(C>G)A", "G(C>G)C", "G(C>G)G", "G(C>G)T",
#               "T(C>G)A", "T(C>G)C", "T(C>G)G", "T(C>G)T",
#               "A(C>T)A", "A(C>T)C", "A(C>T)G", "A(C>T)T",
#               "C(C>T)A", "C(C>T)C", "C(C>T)G", "C(C>T)T",
#               "G(C>T)A", "G(C>T)C", "G(C>T)G", "G(C>T)T",
#               "T(C>T)A", "T(C>T)C", "T(C>T)G", "T(C>T)T",
#               "A(T>A)A", "A(T>A)C", "A(T>A)G", "A(T>A)T",
#               "C(T>A)A", "C(T>A)C", "C(T>A)G", "C(T>A)T",
#               "G(T>A)A", "G(T>A)C", "G(T>A)G", "G(T>A)T",
#               "T(T>A)A", "T(T>A)C", "T(T>A)G", "T(T>A)T",
#               "A(T>C)A", "A(T>C)C", "A(T>C)G", "A(T>C)T",
#               "C(T>C)A", "C(T>C)C", "C(T>C)G", "C(T>C)T",
#               "G(T>C)A", "G(T>C)C", "G(T>C)G", "G(T>C)T",
#               "T(T>C)A", "T(T>C)C", "T(T>C)G", "T(T>C)T",
#               "A(T>G)A", "A(T>G)C", "A(T>G)G", "A(T>G)T",
#               "C(T>G)A", "C(T>G)C", "C(T>G)G", "C(T>G)T",
#               "G(T>G)A", "G(T>G)C", "G(T>G)G", "G(T>G)T",
#               "T(T>G)A", "T(T>G)C", "T(T>G)G", "T(T>G)T")
# # reversed because they go in y axis
# rev_trinuc_96 = rev(trinuc_96)
# 
# maxK = 12
# range_nFact = seq(3, maxK)
# range_k = seq(3, maxK)
# optimal_nFact_k_list = expand.grid(range_nFact, # nFact
#                                    range_k) %>% # K
#   as.matrix %>% t %>% data.frame %>% as.list
# 
# 
# exposures_list = list()
# weights_list = list()
# #reconstructed_list = list()
# 
# exposures_plot_list = list()
# weights_plot_list = list()
# 
# for(optimal_nFact_k in optimal_nFact_k_list){
# 
#   nFact = optimal_nFact_k[1]
#   optimal_k = optimal_nFact_k[2] # the final number of signatures
# 
#   exposures_list[[paste(optimal_nFact_k, collapse = "_")]] = list()
#   weights_list[[paste(optimal_nFact_k, collapse = "_")]] = list()
#   #reconstructed_list[[paste(optimal_nFact_k, collapse = "_")]] = list()
# 
#   for(genome_half in c("1st_half", "2nd_half", "null")){
# 
#     wMatHere = nmfWmatAllByFact_halves[[genome_half]][[sprintf("nFact=%03d", nFact)]]
#     hMatHere = nmfHmatAllByFact_halves[[genome_half]][[sprintf("nFact=%03d", nFact)]] # the NMF w-matrices had been transposed, so rows in the H-matrices and in the W-matrices are the same
#     #reconstructed_list[[paste(optimal_nFact_k, collapse = "_")]][[genome_half]] = list()
#     
#     # ## reconstruct coefficients 
#     # for(actual_nIter in 1:totalNumIters){
#     #   
#     #   actual_nIter = paste0("run", sprintf("%03d", actual_nIter), "_")
#     #   
#     #   exposures_this_actual_nIter = wMatHere %>% 
#     #     data.frame %>% 
#     #     rownames_to_column("nIter_nFact") %>% 
#     #     filter(str_detect(nIter_nFact, actual_nIter)) %>% 
#     #     column_to_rownames("nIter_nFact") %>% 
#     #     as.matrix %>% 
#     #     t
#     #   
#     #   weights_this_actual_nIter = hMatHere %>% 
#     #     data.frame %>% 
#     #     rownames_to_column("nIter_nFact") %>% 
#     #     filter(str_detect(nIter_nFact, actual_nIter)) %>% 
#     #     column_to_rownames("nIter_nFact") %>% 
#     #     as.matrix
#     #   
#     #   reconstructed_list[[paste(optimal_nFact_k, collapse = "_")]][[genome_half]][[actual_nIter]] = exposures_this_actual_nIter %*% weights_this_actual_nIter
#     # }
#     
#     rownames(wMatHere) = gsub("run", "Run ", rownames(wMatHere))
#     rownames(wMatHere) = gsub("fact0", "", rownames(wMatHere))
#     rownames(wMatHere) = gsub("fact", "", rownames(wMatHere))
#     rownames(wMatHere) = gsub("_", "\nFactor ", rownames(wMatHere))
#     rownames(hMatHere) = gsub("run", "Run ", rownames(hMatHere))
#     rownames(hMatHere) = gsub("fact0", "", rownames(hMatHere))
#     rownames(hMatHere) = gsub("fact", "", rownames(hMatHere))
#     rownames(hMatHere) = gsub("_", "\nFactor ", rownames(hMatHere))
# 
#     # need to run the clustering again...
#     set.seed(42 + optimal_k) # this should reproduce the same output as in the "optimizeClustersPar"
#     clustHere = pam(cosineDist(hMatHere), k = optimal_k, diss = T)
#     # medoid: sample in a cluster that is the least dissimilar on average to all the other samples in that cluster
#     signatures = copy(hMatHere[clustHere$medoids,])
#     # sort signatures based on their stability
#     signature_stabilities = data.frame('Signature' = clustHere$medoids, 'Stability' = round(clustHere$silinfo$clus.avg.widths,4)) %>%
#       arrange(Stability)
#     signature_stabilities_sorted = as.character(signature_stabilities$Stability)
#     signatures_sorted = signature_stabilities$Signature
# 
#     # get sample exposures for these cluster medoids (signatures)
#     sample_exposures = wMatHere %>%
#       as_tibble() %>%
#       mutate(signature = rownames(wMatHere)) %>%
#       filter(signature %in% signatures_sorted) %>%
#       arrange(signature) %>% relocate(signature)
# 
#     # and again, [results of pos] + -[results of neg] to have a "-" sign when the larger/non-zero value comes from the neg. submatrix
#     signatures_combined_pos_negcoeff = signatures %>%
#       as_tibble(rownames = NA) %>%
#       rownames_to_column("id") %>%
#       pivot_longer(cols = contains("coeff"),
#                    names_to = "dna_repair_feature_name",
#                    values_to = "Weight") %>%
#       mutate(Weight = ifelse(str_detect(dna_repair_feature_name, "_negcoeff"),
#                              -Weight,
#                              Weight),
#              dna_repair_feature_name = gsub("_...coeff", "", dna_repair_feature_name)) %>%
#       group_by(dna_repair_feature_name, id) %>%
#       summarise(Weight = .Primitive("+")(Weight[1], Weight[2])) %>%
#       ungroup() %>%
#       ## go to original format
#       pivot_wider(names_from = dna_repair_feature_name, values_from = Weight) %>%
#       # arrange iterations (rows)
#       arrange(id) %>%
#       column_to_rownames('id') %>%
#       # arrange dna repair feature_name names (columns)
#       select(coefficient_Resamp_halves[[genome_half]]$feature_name %>%
#                unique) %>%
#       as.matrix()
# 
# 
#     ## parse signature weights
# 
#     # SBS weights in signatures
#     weights = signatures_combined_pos_negcoeff %>%
#       as_tibble() %>%
#       mutate(Signature = rownames(signatures_combined_pos_negcoeff)) %>%
#       left_join(signature_stabilities) %>%
#       pivot_longer(cols = !contains("Signature") & !contains('Stability'), names_to = "feature", values_to = "Weight") %>%
#       mutate(feature_group = ifelse(str_detect(feature, ">"),
#                                     "SBS",
#                                     "Regional\nfeature"),
#              `SBS group` = ifelse(feature_group == "SBS",
#                                   gsub("^.\\(|\\).$", "", feature),
#                                   NA),
#              `SBS group` = factor(`SBS group`,
#                                   levels = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
#                                   ordered = T),
#              SBS96 = ifelse(feature_group == "SBS",
#                             feature,
#                             NA),
#              SBS96 = factor(SBS96,
#                             levels = rev_trinuc_96,
#                             ordered = T),
#              `Regional\nfeature` = ifelse(feature_group == "Regional\nfeature",
#                                           feature,
#                                           NA),
#              # convert signature's stabilities into a factor to be used in y axis, instead of Signature's name which is already indicated in exposures heatmap's y axis
#              Stability = as.character(Stability),
#              Stability = factor(Stability,
#                                 levels = unique(signature_stabilities_sorted))) %>%
#       arrange(SBS96, `Regional\nfeature`, Stability) %>% 
#       mutate("genome_half" = genome_half,
#              "nFact" = nFact,
#              "nK" = optimal_k)
# 
#     weights_list[[paste(optimal_nFact_k, collapse = "_")]][[genome_half]] = weights
# 
# 
#     #### parse signature exposures in samples
#     exposures = sample_exposures %>%
#       pivot_longer(cols = !contains("signature") & !(contains("Stability")), names_to = "Sample", values_to = "Exposure") %>%
#       mutate(Sample = gsub("_nIter.*$", "", Sample)) %>%
#       # add metadata info (e.g. treatments, MSI, HR, smoking...)
#       left_join(samples_info) %>%
#       rename("Signature" = "signature") %>%
#       filter(!is.na(Exposure)) %>%
#       mutate(Signature = factor(Signature, levels = signatures_sorted),
#              "genome_half" = genome_half,
#              "nFact" = nFact,
#              "nK" = optimal_k)
# 
#     exposures_list[[paste(optimal_nFact_k, collapse = "_")]][[genome_half]] = exposures
#   }
# 
#   exposures_plot_list[[paste(optimal_nFact_k, collapse = "_")]] = exposures_list[[paste(optimal_nFact_k, collapse = "_")]] %>% 
#     map(~select(.x, Signature, Sample, genome_half, nFact, nK, Exposure)) %>% 
#     bind_rows() 
#     
#   weights_plot_list[[paste(optimal_nFact_k, collapse = "_")]] = weights_list[[paste(optimal_nFact_k, collapse = "_")]] %>% 
#     map(~select(.x, Signature, feature, genome_half, nFact, nK, Weight)) %>% 
#     bind_rows() 
# }
# 
# 
# 
# ### now compare both halves, should be similar enough, while each very different from the shifted one
# 
# 
# ## EXPOSURES
# 
# exposures_table = bind_rows(exposures_plot_list)
# write_tsv(exposures_table, "heatmap_Ks/exposures_table.tsv")
# 
# ## 1st vs 2nd
# # prepare for plot
# exposures_plot_table = exposures_table %>% 
#   # not using "null" now
#   filter(genome_half != "null") %>% 
#   select(Signature, genome_half, nFact, nK, Exposure) %>% 
#   group_by(Signature, genome_half, nFact, nK) %>% 
#   arrange(Exposure) %>% 
#   mutate(rank = row_number()) %>% 
#   ungroup %>% 
#   pivot_wider(names_from = genome_half, values_from = Exposure) %>% 
#   ##filter out exposures <1
#   #filter_at(vars(contains("_half")), all_vars(.>=1)) %>% 
#   mutate(nFact = factor(gsub("^", "nFact ", nFact),
#                         levels = paste0("nFact ", seq(3,maxK))),
#          nK = factor(gsub("^", "K ", nK),
#                      levels = paste0("K ", seq(maxK,3))))
# 
# # spearman pvals
# exposures_plot_table_spearman = exposures_plot_table %>% 
#   unite('nFact_nK', nFact, nK, sep = "__") %>% 
#   split(as.factor(.$nFact_nK)) %>% 
#   map(~dgof::spearman.test(.x$`1st_half`, .x$`2nd_half`, alternative = "two.sided")) %>% 
#   map(~broom::tidy(.x)) %>% 
#   map(~pull(.x, p.value)) %>% 
#   unlist %>% 
#   data.frame %>% 
#   rownames_to_column("nFact_nK") %>% 
#   separate(nFact_nK, into = c("nFact", "nK"), sep = "__") %>% 
#   rename("K-S p-value" = ".") %>% 
#   mutate(nFact = factor(nFact, levels = paste0("nFact ", seq(3,maxK))),
#          nK = factor(nK, levels = paste0("K ", seq(maxK,3))))
# 
# exposures_plot_table = left_join(exposures_plot_table,
#                                  exposures_plot_table_spearman)
# 
# ggplot(exposures_plot_table,
#        aes(x = `1st_half`,
#            y = `2nd_half`,
#            col = `K-S p-value`)) +
#   geom_point(alpha = 0.5,
#              # shape = 21,
#              # stroke = 0.1,
#              # col = "lightgray",
#              size = 0.5) +
#   # scale_x_log10() +
#   # scale_y_log10() +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   geom_tile() +
#   scale_color_gradientn(colors = c("blue", "lightblue", "red"),
#                         values = c(0, 0.05, 1)) +
#   geom_label(data = exposures_plot_table_spearman,
#              aes(label = round(`K-S p-value`, 2),
#                  x = 0.8, y = 0.2),
#              col = "black",
#              label.padding = unit(0.1, "lines"),
#              label.size = unit(0, "lines"),
#              size = 2) +
#   facet_grid(rows = vars(nK),
#              cols = vars(nFact),
#              scales = "free",
#              space = "free",
#              switch = "both") +
#   xlab("NMF exposures using the ODD chromosomes") +
#   ylab("NMF exposures using the EVEN chromosomes") +
#   theme_classic() +
#   theme(text = element_text(size = 6))
# 
# 
# # prepare for plot
# exposures_plot_table = exposures_table %>% 
#   # not using "null" now
#   filter(genome_half != "null") %>% 
#   select(genome_half, nFact, nK, Exposure) %>% 
#   group_by(genome_half, nFact, nK) %>% 
#   arrange(Exposure) %>% 
#   mutate(rank = row_number()) %>% 
#   ungroup %>% 
#   pivot_wider(names_from = genome_half, values_from = Exposure) %>% 
#   ##filter out exposures <1
#   #filter_at(vars(contains("_half")), all_vars(.>=1)) %>% 
#   mutate(nFact = factor(gsub("^", "nFact ", nFact),
#                         levels = paste0("nFact ", seq(3,maxK))),
#          nK = factor(gsub("^", "K ", nK),
#                      levels = paste0("K ", seq(maxK,3))))
# 
# # spearman pvals
# exposures_plot_table_spearman = exposures_plot_table %>% 
#   unite('nFact_nK', nFact, nK, sep = "__") %>% 
#   split(as.factor(.$nFact_nK)) %>% 
#   map(~dgof::spearman.test(.x$`1st_half`, .x$`2nd_half`, alternative = "two.sided")) %>% 
#   map(~broom::tidy(.x)) %>% 
#   map(~pull(.x, p.value)) %>% 
#   unlist %>% 
#   data.frame %>% 
#   rownames_to_column("nFact_nK") %>% 
#   separate(nFact_nK, into = c("nFact", "nK"), sep = "__") %>% 
#   rename("K-S p-value" = ".") %>% 
#   mutate(nFact = factor(nFact, levels = paste0("nFact ", seq(3,maxK))),
#          nK = factor(nK, levels = paste0("K ", seq(maxK,3))))
# 
# exposures_plot_table = left_join(exposures_plot_table,
#                                  exposures_plot_table_spearman)
# 
# ggplot(exposures_plot_table,
#        aes(x = `1st_half`,
#            y = `2nd_half`,
#            col = `K-S p-value`)) +
#   geom_point(alpha = 0.5,
#              # shape = 21,
#              # stroke = 0.1,
#              # col = "lightgray",
#              size = 0.5) +
#   # scale_x_log10() +
#   # scale_y_log10() +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   geom_tile() +
#   scale_color_gradientn(colors = c("blue", "lightblue", "red"),
#                         values = c(0, 0.05, 1)) +
#   geom_label(data = exposures_plot_table_spearman,
#              aes(label = round(`K-S p-value`, 2),
#                  x = 0.8, y = 0.2),
#              col = "black",
#              label.padding = unit(0.1, "lines"),
#              label.size = unit(0, "lines"),
#              size = 2) +
#   facet_grid(rows = vars(nK),
#              cols = vars(nFact),
#              scales = "free",
#              space = "free",
#              switch = "both") +
#   xlab("NMF exposures using the ODD chromosomes") +
#   ylab("NMF exposures using the EVEN chromosomes") +
#   theme_classic() +
#   theme(text = element_text(size = 6))
# 
# 
# ## 1st vs. null
# # prepare for plot
# exposures_plot_table = exposures_table %>% 
#   # not using "2nd" now
#   filter(genome_half != "2nd_half") %>% 
#   select(genome_half, nFact, nK, Exposure) %>% 
#   group_by(genome_half, nFact, nK) %>% 
#   arrange(Exposure) %>% 
#   mutate(rank = row_number()) %>% 
#   ungroup %>% 
#   pivot_wider(names_from = genome_half, values_from = Exposure) %>% 
#   mutate(nFact = factor(gsub("^", "nFact ", nFact),
#                         levels = paste0("nFact ", seq(3,maxK))),
#          nK = factor(gsub("^", "K ", nK),
#                      levels = paste0("K ", seq(maxK,3))))
# 
# # spearman pvals
# exposures_plot_table_spearman = exposures_plot_table %>% 
#   unite('nFact_nK', nFact, nK, sep = "__") %>% 
#   split(as.factor(.$nFact_nK)) %>% 
#   map(~dgof::spearman.test(.x$`1st_half`, .x$`null`, alternative = "two.sided")) %>% 
#   map(~broom::tidy(.x)) %>% 
#   map(~pull(.x, p.value)) %>% 
#   unlist %>% 
#   data.frame %>% 
#   rownames_to_column("nFact_nK") %>% 
#   separate(nFact_nK, into = c("nFact", "nK"), sep = "__") %>% 
#   rename("K-S p-value" = ".") %>% 
#   mutate(nFact = factor(nFact, levels = paste0("nFact ", seq(3,maxK))),
#          nK = factor(nK, levels = paste0("K ", seq(maxK,3))))
# 
# exposures_plot_table = left_join(exposures_plot_table,
#                                  exposures_plot_table_spearman)
# 
# ggplot(exposures_plot_table,
#        aes(x = `1st_half`,
#            y = `null`,
#            col = `K-S p-value`)) +
#   geom_point(alpha = 0.5,
#              size = 0.5) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   geom_tile() +
#   scale_color_gradientn(colors = c("blue", "lightblue", "red"),
#                         values = c(0, 0.05, 1)) +
#   geom_label(data = exposures_plot_table_spearman,
#              aes(label = round(`K-S p-value`, 2),
#                  x = 0.8, y = 0.2),
#              col = "black",
#              label.padding = unit(0.1, "lines"),
#              label.size = unit(0, "lines"),
#              size = 2) +
#   facet_grid(rows = vars(nK),
#              cols = vars(nFact),
#              scales = "free",
#              space = "free",
#              switch = "both") +
#   xlab("NMF exposures using the ODD chromosomes") +
#   ylab("NMF exposures using the NULL data") +
#   theme_classic() +
#   theme(text = element_text(size = 6))
# 
# 
# 
# ## WEIGHTS
# 
# weights_plot_table = bind_rows(weights_plot_list)
# write_tsv(weights_plot_table, "heatmap_Ks/weights_plot_table.tsv")
# 
# ggplot(weights_plot_table)
