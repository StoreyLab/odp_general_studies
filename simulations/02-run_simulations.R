#! /usr/bin/env Rscript
source("./02-main_simulations.R")
library(reshape2)
library(edge)
library(qvalue)
library(dplyr)
library(Biobase)
library(splines)
library(qvalue)
args <- commandArgs(trailingOnly = TRUE)
experiment <- args[1]
num.curves <- as.numeric(args[2])
replicates <- as.numeric(args[3])
seed <- as.numeric(args[4])
infolder <- args[5]
outfolder <- args[6]
print(seed)
print(args)
ts <- run_simulation(experiment = experiment,
                             num.curves = num.curves, infolder = infolder, outfolder = outfolder,
                             replicates = replicates,
                             seed = seed)

save(ts, file = paste0(outfolder, "rdafiles/", experiment, "_", num.curves, "_", replicates, ".rda"))
