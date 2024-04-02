library(tidyverse)

permuted_coeff_null = read_tsv("permuted_coefficients_100iters_null_ORIG.tsv")

permuted_coeff_null_feat_rm = permuted_coeff_null %>% 
  dplyr::filter(!feature_name %in% c("MBD4","SETD2","53BP1","KAT2A","CREBBP"))

write_tsv(permuted_coeff_null_feat_rm,
          "permuted_coefficients_100iters_null.tsv")
