source("helper.R")
library(tidyverse)
library(edge)
library(biobroom)

###############################################################################
# Wrapper of estFDR for each simulation setting
# infile: string, location of simulation file
###############################################################################
summarise_simulations <- function(infile) {
  load(infile)
  ts %>%
    dplyr::group_by(nc, variable, replicate, snr_sampling, null_curve) %>%
    dplyr::do(estFDR(oracle = .$oracle, pvals = .$value))
}

###############################################################################
# Function to summarise simulations (estimated FDR, total sig genes)
# infile: string, location of simulation file
###############################################################################
run_summarise <- function(infolder, outfolder) {
  infiles <- list.files(infolder, full.names = TRUE)
  fname <- list.files(infolder)
  dir.create(outfolder, showWarnings = FALSE)
  cmds <- paste("sbatch -p storey",  "~/job_timecourse_summarise.sh", infiles, fname, outfolder)
  for (cmd in cmds) {
    cat(cmd, "\n")
    system(cmd)
  }
}

###############################################################################
# Execute time course simulations
# experiment: experiment name
# num.curves: number of unique time profiles
# replicates: number of simulation replicates
# infolder: location of simulation file
# outfolder: location of output folder
# seed: seed value
###############################################################################
run_simulation <- function(experiment, num.curves, replicates, infolder, outfolder, seed) {
  print(seed)
  set.seed(seed)
  cutoff <- ifelse(experiment == "endotoxin",  .05,  .1)
  obj <- readRDS(paste0(infolder, experiment, "_lrt.rds"))
  snr <- get_SNR(obj, fdr_cutoff = cutoff)
  fit <- edge::fit_models(obj, stat.type = "lrt")
  q.thres <- ifelse(experiment == "endotoxin",  .45,  ifelse(experiment == "smoker", .8, ifelse(experiment == "kidney", .35, .86)))
  # Model fits for null/alt
  fit.full <- fit@fit.full
  fit.null <- fit@fit.null
  dims_eo <- dim(obj)
  # Estimated pi0 and expected number of alternatives
  pi0 <- obj@qvalueObj$pi0
  totSigGenes <- round((1 - pi0) * dims_eo[1] / num.curves) * num.curves
  num_nulls <- dims_eo[1] - totSigGenes
  simulateData <- function(num.curves, experiment) {
    select_curves <- sample(snr$alt_id, size = num.curves, replace = FALSE)
    null_ind <- sample(snr$null_id, size = num_nulls, replace = TRUE)
    sig_genes <- fit.full[sample(select_curves, replace = TRUE, size = totSigGenes),]
    null_genes <- fit.null[null_ind,]
    edata <- rbind(sig_genes, null_genes)
    snr <- get_SNR(obj, fdr_cutoff = 1)
    noise_sd <- apply(sig_genes, 1, sd) / sqrt(quantile(snr$snr, q.thres))
    noise_sd <- c(noise_sd, sample(noise_sd, replace = TRUE, size = num_nulls))
    noise <- t(mapply(noise_sd, FUN = function(x) rnorm(ncol(edata), mean = 0, sd = x)))
    return(edata + noise)
  }
  sim_data <- simulateData(num.curves, experiment = experiment)
  status <- c(rep(1, totSigGenes), rep(0, num_nulls))
  rownames(sim_data) <- rownames(obj)
  colnames(sim_data) <- colnames(obj)
  exprs(obj) <- as.matrix(sim_data)
  if (experiment == "endotoxin") {
    individual(obj) <- as.factor(NULL)
  }
  
  rownames(obj) <- 1:nrow(obj)
  colnames(obj) <- 1:ncol(obj)
  odp <- augment(odp(obj, bs.its = 500, n.mods = 800, seed = seed + 1)@qvalueObj,
          Method = "ODP",
          seed = seed + 1,
          study = experiment,
          nc = num.curves, 
          replicate = replicates, oracle = status)
  
  lrt <- augment(lrt(obj, seed = seed + 2)@qvalueObj,
                 Method = "LRT",
                 seed = seed + 2,
                 study = experiment,
                 nc = num.curves, 
                 replicate = replicates, oracle = status)
  
  mod_lrt <- augment(lrt(obj, mod.F = TRUE, seed = seed + 3)@qvalueObj,
                 Method = "mLRT",
                 seed = seed + 3,
                 study = experiment,
                 nc = num.curves, 
                 replicate = replicates, oracle = status)
  
  blrt <- augment(lrt(obj, nullDistn = "bootstrap", bs.its = 500, seed = seed + 4)@qvalueObj,
                 Method = "bLRT",
                 seed = seed + 4,
                 study = experiment,
                 nc = num.curves, 
                 replicate = replicates, oracle = status)
  
  bmod_lrt <- augment(lrt(obj, nullDistn = "bootstrap", bs.its = 500, mod.F = TRUE, seed = seed + 5)@qvalueObj,
                     Method = "bmLRT",
                     seed = seed + 5,
                     study = experiment,
                     nc = num.curves, 
                     replicate = replicates, oracle = status)
  
  bind_rows(odp, lrt, mod_lrt, blrt, bmod_lrt)
}



###############################################################################
# Run time course simulations
# outfolder: location to put results
# replicates: total number of replicates
# seed: seed of simulation
# infolder: string, location of simulation file
###############################################################################
timecourse_run <- function(outfolder = "../data/raw_data/", infolder = "../analysis/data/", replicates = 1:50, seed = 123) {
  # initializations for simulation
  dir.create(outfolder, showWarnings = FALSE)
  curves <- c(5, 10, 50, 100, 200)
  exp <- c("endotoxin", "kidney", "dose", "smoker")
  
  designs <- expand.grid(experiment = exp,
                      num.curves = curves,
                      replicates = replicates)

  set.seed(seed)
  designs$seed <- sample(.Machine$integer.max, size = nrow(designs))
  cmds <- paste("sbatch -p storey",  "./job_simulations.sh", designs$experiment,
                designs$num.curves, designs$replicates, designs$seed, infolder, outfolder)
  print(cmds)
  for (cmd in cmds) {
    system(cmd)
  }
}
