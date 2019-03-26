estFDR <- function(oracle, qvals, fdr = seq(0, .1, 0.001)) {
  est_fdr <- tot_sigGenes <- rep(0, length(fdr))
  for (i in 1:length(fdr)) {
    est <- oracle[qvals < fdr[i]]
    tot_sigGenes[i] <- length(est)
    FP = sum(!est)
    est_fdr[i] <- FP / pmax(length(est), 1)
  }
  return(data.frame(true_fdr = fdr,
                    tot_genes = tot_sigGenes,
                    est_fdr = est_fdr))
}

summarise_simulations <- function(infile) {
  load(infile)
  ts %>%
    dplyr::group_by(nc, replicate, study, Method) %>%
    dplyr::do(estFDR(oracle = .$oracle, qvals = .$q.value))
}

library(tidyverse)
ll <- list.files("../data/raw_data/rdafiles/")
ts  <- NULL
for (i in 1:length(ll)) {
  print(i)
  ts <- rbind(ts,summarise_simulations(infile = paste0("../data/raw_data/rdafiles/", ll[i])))
}

saveRDS(ts, "../data/simulation_results.rds")
