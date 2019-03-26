###############################################################################
# Analyze endotoxin data
# seed: seed of simulation
# infile: string, location of raw data
###############################################################################
analyze_endotoxin <- function(infile = "../data/raw_data/", outfile = "../data/raw_data/", seed = 12345) {
  require(SummarizedExperiment)
  require(splines)
  edata <- read.table(paste0(infile, "CR001_expr.txt"),
                      header = FALSE)
  pdata <- read.table(paste0(infile, "CR001_cov.txt"),
                      header = FALSE,
                      row.names = 1)
  gnames = read.table(paste0(infile, "CR001_genes.txt"))

  # transpose phenotype data
  pdata <- as_data_frame(t(pdata))
  colnames(pdata) <- c("class", "individual", "time")
  
  # Change column type
  pdata$time <- as.numeric(pdata$time)
  pdata$individual <- as.factor(pdata$individual)
  pdata$class <- as.factor(pdata$class) 
  
  # transform expression data
  edata = log2(as.matrix(edata) + 10)
  
  # create models
  null_model <- ~ns(time, df = 4, intercept = FALSE) 
  full_model <- ~ns(time, df = 4, intercept = FALSE) * class 

  colnames(edata) <- rownames(pdata)
  de_obj <- edge::build_models(data = edata,
                               cov = pdata,
                               full.model = full_model,
                               null.model = null_model,
                               ind = pdata$individual)
  obj <- edge::lrt(de_obj, nullDistn = "bootstrap", bs.its = 500, seed = seed)
  saveRDS(obj, file = paste0(outfile, "endotoxin_lrt_boot.rds"))
  obj <- edge::lrt(de_obj, seed = seed + 1)
  saveRDS(obj, file = paste0(outfile, "endotoxin_lrt.rds"))
  obj <- edge::lrt(de_obj, nullDistn = "bootstrap", bs.its = 500, seed = seed+2)
  saveRDS(obj, file = paste0(outfile, "endotoxin_mlrt_boot.rds"))
  obj <- edge::lrt(de_obj, mod.F = TRUE, seed = seed+3)
  saveRDS(obj, file = paste0(outfile, "endotoxin_mlrt.rds"))
  obj <- edge::odp(de_obj,  n.mods = 800, bs.its = 500, seed = seed+4)
  saveRDS(obj, file = paste0(outfile, "endotoxin_odp.rds"))
  obj
}


###############################################################################
# Analyze kidney data
# seed: seed of simulation
# infile: string, location of raw data
###############################################################################
analyze_kidney <- function(infile = "../data/raw_data/", outfile = "../data/raw_data/", seed = 12345) {
  require(tidyverse)
  require(splines)
  edata <- read_delim(paste0(infile, "kidney_expr.txt"), delim = "\t", col_names = FALSE)
  pdata <-  t(read_delim(paste0(infile, "kidney_cov.txt"), delim = "\t", col_names = FALSE))
  filter.probes <- readRDS(paste0(outfile, "filter_probes.rds"))
  colnames(pdata) <- c("tissue", "age","pathology", "creatinine", "sex","gfr")
  pdata = type_convert(as_tibble(pdata[-1,]))
  
  # Only look at cortex samples 
  edata <- log2(edata[filter.probes, pdata$tissue == "c"] + 10)
  pdata <- pdata[pdata$tissue == "c",]
  pdata$sex <- as.factor(pdata$sex)
  colnames(edata) <- rownames(pdata) <- NULL
  
  # create models
  full_model <- ~sex + ns(age, df = 4, intercept = FALSE)
  null_model <- ~sex
  
  de_obj <- edge::build_models(data = as.matrix(edata),
                               cov = pdata,
                               full.model = full_model,
                               null.model = null_model)

  obj <- edge::lrt(de_obj, nullDistn = "bootstrap", bs.its = 500, seed = seed)
  saveRDS(obj, file = paste0(outfile, "kidney_lrt_boot.rds"))
  obj <- edge::lrt(de_obj, seed = seed + 1)
  saveRDS(obj, file = paste0(outfile, "kidney_lrt.rds"))
  obj <- edge::lrt(de_obj, nullDistn = "bootstrap", bs.its = 500, seed = seed+2)
  saveRDS(obj, file = paste0(outfile, "kidney_mlrt_boot.rds"))
  obj <- edge::lrt(de_obj, mod.F = TRUE, seed = seed+3)
  saveRDS(obj, file = paste0(outfile, "kidney_mlrt.rds"))
  obj <- edge::odp(de_obj, n.mods = 800, bs.its = 500, seed = seed+4)
  saveRDS(obj, file = paste0(outfile, "kidney_odp.rds"))
  obj
}

###############################################################################
# Analyze dose-response data
# seed: seed of simulation
# infile: string, location of raw data
###############################################################################
analyze_dose <- function(infile = "../data/raw_data/", outfile = "../data/raw_data/", seed = 12345) {
  require(SummarizedExperiment)
  require(GEOquery)
  require(stringr)
  require(tidyverse)
  
  gds <- getGEO("GSE4668") 
  eset <- gds$GSE4668_series_matrix.txt.gz

  pdata <- as.data.frame(eset@phenoData@data) %>%
    select(title, geo_accession) %>%
    mutate(title = as.character(title))
  pdata <- pdata[1:25,]
  pdata$dose <- as.numeric(unlist(lapply(str_split(pdata$title, " "), FUN = function(x) x[2]))[1:25])
  expression <- log2(eset@assayData$exprs[, 1:25] + 10)

  # Model
  fmod = ~1 + ns(dose, df = 2)
  nmod = ~1
  rownames(expression) <-  colnames(expression) <- rownames(pdata) <- NULL

  de_obj <- edge::build_models(data = as.matrix(expression),
                               cov = pdata,
                               full.model = fmod,
                               null.model = nmod)

  obj <- edge::lrt(de_obj, nullDistn = "bootstrap", bs.its = 500, seed = seed)
  saveRDS(obj, file = paste0(outfile, "dose_lrt_boot.rds"))
  obj <- edge::lrt(de_obj,  seed = seed + 1)
  saveRDS(obj, file = paste0(outfile, "dose_lrt.rds"))
  obj <- edge::lrt(de_obj, nullDistn = "bootstrap", mod.F= TRUE,bs.its = 500, seed = seed+2)
  saveRDS(obj, file = paste0(outfile, "dose_mlrt_boot.rds"))
  obj <- edge::lrt(de_obj,  mod.F = TRUE, seed = seed+3)
  saveRDS(obj, file = paste0(outfile, "dose_mlrt.rds"))
  obj <- edge::odp(de_obj,  n.mods = 800, bs.its = 500, seed = seed+4)
  saveRDS(obj, file = paste0(outfile, "dose_odp.rds"))
  obj
}

###############################################################################
# Analyze RNA-Seq static smoker data
# seed: seed of simulation
# infile: string, location of raw data
###############################################################################
analyze_smoker <- function(infile = "../data/raw_data/", outfile = "../data/raw_data/", seed = 12345) {
  require(stringr)
  require(splines)
  require(edge)
  require(tidyverse)
  require(ExpressionAtlas)
  experimentSummary <- getAtlasExperiment(experimentAccession  = "E-GEOD-47718")
  rse_gene = experimentSummary$rnaseq
  edata <- assay(rse_gene)
  edata = edata[rowSums(edata) >= 10, ]
  
  pdata <- colData(rse_gene)
  
  # Get covariates
  pdata$condition = pdata$clinical_information
  
  pdata$technical_replicate_group[pdata$technical_replicate_group[1] == pdata$technical_replicate_group] <- rep(paste0("group", 17:1), each = 1)
  colnames(edata) <- pdata$technical_replicate_group
  edata = reshape2::melt(edata) %>% group_by(Var1, Var2) %>% summarise(total = sum(value)) %>% spread(Var2, total)
  edata = edata[, -1]
  pdata = unique(pdata)
  
  # create models
  null_model <- ~1
  full_model <- ~factor(condition)
  pdata <- as.data.frame(pdata) %>% select(condition)
  
  # mean-variance relationship
  vm = limma::voom(edata, design = model.matrix(full_model, pdata))
  colnames(vm$E) <- rownames(pdata)

  de_obj <- edge::build_models(data = vm$E,
                               cov = pdata,
                               full.model = full_model,
                               null.model = null_model)
  
  weights = vm$w
  obj <- edge::lrt(de_obj, weights = weights, nullDistn = "bootstrap", bs.its = 5000, seed = seed)
  saveRDS(obj, file = paste0(outfile, "smoker_lrt_boot.rds"))
  obj <- edge::lrt(de_obj, weights = weights, seed = seed + 1)
  saveRDS(obj, file = paste0(outfile, "smoker_lrt.rds"))
  obj <- edge::lrt(de_obj, weights = weights, nullDistn = "bootstrap", mod.F = TRUE, bs.its = 5000, seed = seed+2)
  saveRDS(obj, file = paste0(outfile, "smoker_mlrt_boot.rds"))
  obj <- edge::lrt(de_obj, weights = weights, mod.F = TRUE, seed = seed+3)
  saveRDS(obj, file = paste0(outfile, "smoker_mlrt.rds"))
  obj <- edge::odp(de_obj, weights = weights, n.mods = 800, bs.its = 5000, seed = seed+4)
  saveRDS(obj, file = paste0(outfile, "smoker_odp.rds"))
  obj
}

library(reshape2)
library(edge)
library(qvalue)
library(dplyr)
library(Biobase)
library(splines)
library(qvalue)
args <- commandArgs(trailingOnly = TRUE)
experiment <- args[1]
print(args)
if (experiment == "dose") analyze_dose()
if (experiment == "endotoxin") analyze_endotoxin()
if (experiment == "kidney") analyze_kidney()
if (experiment == "smoker") analyze_smoker()
