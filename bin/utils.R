library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("between", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("reduce", "IRanges")
conflict_prefer("desc", "dplyr")



# Overlaps wrapper; each "ranges_[A,B]" is a tracks file stratified by chromosome as a list of granges
findOverlaps_wrapper <- function(ranges_A, ranges_B) {

  ## get queryHits & subjectHits
  findOverlaps_df = findOverlaps(query = ranges_A,
                                 subject = ranges_B) %>% 
    data.frame

  ## convert ranges_A and _B to df; ONLY the ranges that overlap, i.e. queryHits for A and subjectHits for B
  ranges_A_df = data.frame(ranges_A[findOverlaps_df$queryHits])
  ranges_B_df = data.frame(ranges_B[findOverlaps_df$subjectHits])

  ## append ranges_A_df and _B_df $start and $end to queryHits table
  output = findOverlaps_df %>% 
    select(queryHits) %>% 
    mutate(ranges_A_start = as.numeric(ranges_A_df$start),
           ranges_A_end = as.numeric(ranges_A_df$end),
           ranges_B_start = as.numeric(ranges_B_df$start),
           ranges_B_end = as.numeric(ranges_B_df$end)) %>% 
    # also append the chr and scores binary levels
    bind_cols(select(ranges_A_df,
                     -c(start, end, width, strand))) %>% 
    relocate(seqnames, .after = queryHits) %>% 
    ## finish getting overlapping ranges correct 'start' and 'end' values
    rowwise() %>% 
    mutate(# highest 'start' value
           start = max(ranges_A_start, ranges_B_start),
           # lowest 'end' value
           end = min(ranges_A_end, ranges_B_end)) %>%
    relocate(start, .after = seqnames) %>% relocate(end, .after = start) %>% 
    select(-c(queryHits, contains("ranges_"))) %>%
    ungroup() %>% 
    arrange(start)
}



###########################################################################################################################

### match trinucleotide(32) proportions between genomic coordinates
# takes trinucleotide(32) count dataframe with 2 rows (1 per bin (high/low)) × 32 trinuc type columns
# removes trinucs(32) counts, so that trinucs(32) proportions are "matched" across all bins
# modified from data/users/malvarez/projects/RepDefSig/resources/trinuc_matching/marina/sampler_fran.R

# normalize each row to sum to 1
rowNorm = function(m){
  m / rowSums(m)
}

# euclidean distance
euclidean = function(a, b){
  sqrt(sum((a-b)^2))
} 

trinuc_matching = function(trinuc_matching_input, 
                           stoppingCriterion = 0.01, # desired Euclidean score (max. overall distance between any bin's trinuc frequencies and all-bin-average trinuc frequencies)
                           maxIter = 1e6, # to prevent endless loops
                           maxTime = 1, # max time
                           unitsTime = "hours", # time units, if >24h should specify "days"
                           max_fraction_removed_trinucs = 0.2, # don't allow to remove more total trinucleotide counts than this fraction of the total original trinucleotide counts
                           fast_progress_lfc_cutoff = -0.00001, # minimum degree of Euclidean score decrease (LFC; e.g. log2(0.0615/0.062)) allowed for the mean of progress_its
                           progress_its = 1000, # n last iterations (that reached a mineuclidean_score) used to calculate the progress
                           seed = 1){
  ## initialize constants/variables
  set.seed(seed)
  euclidean_score = sqrt(2) # we want to make this decrease until 'stoppingCriterion'
  mineuclidean_score = euclidean_score
  progress = c() # to keep track of how the euclidean_score decrease is progressing (and early-stop if it shows signs of stagnation)
  counts = trinuc_matching_input
  min_counts = min(counts)
  counts_mineuclidean = counts
  removed_trinucs = 0
  rowSums_counts = rowSums(counts)
  total_orig_trinucs = sum(rowSums_counts)
  max_removed_trinucs = total_orig_trinucs * max_fraction_removed_trinucs
  iter = 0
  
  # each loop iteration should take not more than 5 minutes. Adjust this value accordingly.
  if(unitsTime == "hours"){
    expected_max_duration_it = 5/60 # 5 minutes in hours
  } else if(unitsTime == "days"){
    expected_max_duration_it = 5/60/24 # 5 minutes in days
  }
  
  total_analysis_time = 0 # analysis will stop if this value reaches maxTime 
  last_checkpoint_time = Sys.time()
  
  while ( TRUE ) {
    
    #### adjust for slurm job suspension events
    time_since_last_checkpoint = as.numeric(difftime(Sys.time(), # current time
                                                     last_checkpoint_time, units = unitsTime))
    
    if (time_since_last_checkpoint > expected_max_duration_it) {
      # Detected a suspension event (more than 5 minutes passed since previous iteration)
      total_analysis_time = total_analysis_time + expected_max_duration_it
    } else {
      total_analysis_time = total_analysis_time + time_since_last_checkpoint
    }
    
    ## check whether we have reached max. time
    if (total_analysis_time >= maxTime){
      counts = counts_mineuclidean
      cat( sprintf("Stopping optimization - %s %s have passed - Returning min Euclidean score results (%f)\n", total_analysis_time, unitsTime, mineuclidean_score) )
      break
    }
    
    # Update last_checkpoint_time
    last_checkpoint_time = Sys.time()
    
    
    ####
    ## continue with iteration
    
    iter = iter + 1
    
    ## check whether we have reached max. num. iterations
    if (iter == maxIter) {
      counts = counts_mineuclidean
      cat( sprintf("Stopping optimization - maximum number of iterations reached - Returning min Euclidean score results (%f)\n", mineuclidean_score) )
      break
    }  
    
    ## check whether we have already removed too many trinucs
    if(removed_trinucs > max_removed_trinucs){
      counts = counts_mineuclidean
      cat(sprintf("Stopping optimization at iter %i/%i - %.02f%% of the %f original trinucs have already been removed -- Exiting and returning min Euclidean score results (%f)\nAnalysis terminated after %s %s\n", iter, maxIter, removed_trinucs/total_orig_trinucs*100, total_orig_trinucs, mineuclidean_score, total_analysis_time, unitsTime))
      break
    }

    
    ###################################################
    ## get offender trinuc and bin
    
    ## identify the trinuc (column) whose difference in relative frequency between the 2 bins is the highest: offender
    trinuc_proportions = rowNorm(counts)
    
    offender_trinuc_proportion_diffs = diff(as.matrix(trinuc_proportions)) %>% 
      data.frame() %>% 
      rownames_to_column("bin")
    
    # which is the reference bin level, i.e. the one from which we subtract the trinuc proportions from the other level
    bin_ref = offender_trinuc_proportion_diffs$bin
    
    offender_trinuc = offender_trinuc_proportion_diffs %>% 
      select(-bin) %>% 
      pivot_longer(cols = everything(), names_to = "trinuc32", values_to = "diff_proportion") %>% 
      slice_max(abs(diff_proportion))
    
    # if there is a tie between >=2 trinucs, randomly select one
    if(length(offender_trinuc) >= 2){ offender_trinuc = offender_trinuc %>% slice_sample(n = 1) }
    
    # this diff_proportion will be used to subsample the offender column counts at the offender bin
    offender_trinuc_diff_proportion = offender_trinuc %>% 
      pull(diff_proportion)
    
    offender_trinuc = offender_trinuc %>% 
      pull(trinuc32)
      
    ## the offender bin is the row with HIGHER PRPORTION of the current offender column (even if it has lower absolute counts), to be subsampled
    if(offender_trinuc_diff_proportion > 0){
      
      # diff_proportion positive, so the bin_ref has the highest relative proportion for this offender column, so counts will be deducted from this bin
      offender_bin_trinuc = counts %>% 
        select(all_of(offender_trinuc)) %>% 
        rownames_to_column("bin") %>% 
        filter(bin == bin_ref)
      
    } else {
      
      # diff_proportion negative, counts deducted from the non-ref. bin
      offender_bin_trinuc = counts %>% 
        select(all_of(offender_trinuc)) %>% 
        rownames_to_column("bin") %>% 
        filter(bin != bin_ref)
    }
    
    offender_bin = offender_bin_trinuc %>% 
      pull(bin)
    
    offender_bin_counts = offender_bin_trinuc %>% 
      pull(!!sym(offender_trinuc))
    
    ###########


    # store euclidean_score from previous round
    previous_euclidean_score = euclidean_score
    
    ## calculate euclidean score (how different (i.e. euclidean distance) are the trinuc freqs. between the 2 bins)
    euclidean_score = euclidean(trinuc_proportions[1,], trinuc_proportions[2,])
    
    # calculate ratio of change in this iteration
    euclidean_change_ratio = euclidean_score / previous_euclidean_score
    
    # if euclidean_score is new minimum:
    if(euclidean_score < mineuclidean_score){
      
      # 1- store counts table
      mineuclidean_score = euclidean_score
      counts_mineuclidean = counts
      
      # 2- calculate amount of change between previous and current iteration (since it reached a mineuclidean_score)...
      LFC_euc_score = log2(euclidean_change_ratio)
      
      # ...and append it to the 'amount_change's of the last 'progress_its' iterations (that reached a mineuclidean_score)
      progress = c(tail(progress, progress_its-1),
                   LFC_euc_score)
      mean_progress_lfc = mean(progress)
    }
    
    # did we reduce the euclidean_score enough?
    if (euclidean_score <= stoppingCriterion) {
      cat(sprintf("Successfully completed optimization: Euclidean score (%f) lower than %f - Returning current results\n", euclidean_score, stoppingCriterion) )
      break
    } else if (mean_progress_lfc>fast_progress_lfc_cutoff  &  length(progress)==progress_its) {
      # if not, is the euclidean_score stagnated? (i.e. very small/slow progress)
      cat(sprintf("Cannot continue optimization at iter %i/%i - progress (mean LFC) of last %i iterations (that decreased the min. Euclidean score) is too slow: %f > %f (fast_progress_lfc_cutoff) -- Exiting and returning min Euclidean score results (%f)\nAnalysis terminated after %s %s\n", iter, maxIter, progress_its, mean_progress_lfc, fast_progress_lfc_cutoff, mineuclidean_score, total_analysis_time, unitsTime))
      break
    }
    

    ###############################
    ### subtract counts from the responsible trinuc (column) of the offender bin (row) so that the offender column frequency of the latter get closer to the offender column frequency of the non-offender bin
    
    subtractThis = round(offender_bin_counts * abs(offender_trinuc_diff_proportion))

    # total counts that will remain at the offender_name row × correctableCol trinuc intersection after the removal of 'subtractThis' counts
    remaining_offender_bin_counts = offender_bin_counts - subtractThis
    
    # make sure that there were not more counts removed than the total available
    if(remaining_offender_bin_counts < 0){
      
      # if so, let there be min_counts (lowest counts bin×trinuc combination) left at that row×col intersection
      remaining_offender_bin_counts = min_counts
        
      # update subtractThis just to keep track of the actual # of counts removed
      subtractThis = min_counts - offender_bin_counts
    }
    
    # keep track of # of counts removed
    removed_trinucs = removed_trinucs + abs(subtractThis)
    
    ### apply this change to the actual counts table
    counts = counts %>% 
      mutate(!!sym(offender_trinuc) := replace(!!sym(offender_trinuc), row_number() == which(rownames(counts) == offender_bin), 
                                               remaining_offender_bin_counts))
    
    ## output log every 50th iter
    if(iter %% 50 == 0){
        
      # print log message
      cat( sprintf("Iteration %i/%i:\n\tSubtracted %i '%s's at bin '%s'\n\t%.02f%% of the original trinucleotides have been removed\n\tEuclidean score: %f\n\tMean progress LFC: %f for the last %i iterations that updated the min. Euclidean score\n\t%s %s have passed\n\n", iter, maxIter, abs(as.numeric(subtractThis)), offender_trinuc, offender_bin, removed_trinucs/total_orig_trinucs*100, euclidean_score, mean_progress_lfc, progress_its, total_analysis_time, unitsTime))
    }
    
    ## save checkpoint every 10K iteration
    if(iter %% 1e4 == 0){
      
      matched_tracks_checkpoint = counts_mineuclidean %>%
        # put back the bin names
        rownames_to_column("bin")
      
      saveRDS(matched_tracks_checkpoint, paste0(as.character(format(iter, scientific = F)), 
                                                "th_iter__euclidean_", 
                                                round(mineuclidean_score,2), 
                                                "__matched_tracks_checkpoint.Rda"))
    }
  } ## keep iterating...
  
  ## return final counts_mineuclidean tables, and its mineuclidean_score
  return(list(counts_mineuclidean %>%
                # put back the bin names
                rownames_to_column("bin"),
              mineuclidean_score))
}
