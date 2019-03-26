library(biobroom)
library(edge)
require(ExpressionAtlas)
library(tidyverse)
require(SummarizedExperiment)
require(GEOquery)
require(stringr)

#### Combine kidney ####
outfile <- "../analysis/data/"
filter.probes <- readRDS(paste0(outfile, "filter_probes.rds"))
rnames <- read_delim("../analysis/data/raw_data/kidney_genes.txt", delim = " ", col_names = FALSE)
rnames <- rnames$X1[filter.probes]
o <- readRDS("../analysis/data/kidney_odp.rds")
o <- augment(o@qvalueObj, Method = "ODP", names = rnames, study = "Kidney")
l <- readRDS("../analysis/data/raw_data/kidney_lrt.rds")
l <- augment(l@qvalueObj, Method = "F-test", names = rnames, study = "Kidney")
lm <- readRDS("../analysis/data/raw_data/kidney_mlrt.rds")
lm <- augment(lm@qvalueObj, Method = "Moderated\nF-test", names = rnames, study = "Kidney")
l2 <- readRDS("../analysis/data/raw_data/kidney_lrt_boot.rds")
l2 <- augment(l2@qvalueObj, Method = "boot F-test", names = rnames, study = "Kidney")
lm2 <- readRDS("../analysis/data/raw_data/kidney_mlrt_boot.rds")
lm2 <- augment(lm2@qvalueObj, Method = "boot Moderated\nF-test", names = rnames, study = "Kidney")

df_kidney <- rbind(o, l, lm, l2, lm2)
saveRDS(df_kidney, file = "../analysis/data/kidney_comp.rds")


##### Smoker dataset
experimentSummary <- getAtlasExperiment(experimentAccession  = "E-GEOD-47718")
rse_gene = experimentSummary$rnaseq
edata <- assay(rse_gene)
edata = edata[rowSums(edata) >= 10, ]
rnames <- rownames(edata)
o <- readRDS("../analysis/data/raw_data/smoker_odp.rds")
o <- augment(o@qvalueObj, Method = "ODP", names = rnames, study = "Smoker")
l <- readRDS("../analysis/data/raw_data/smoker_lrt.rds")
l <- augment(l@qvalueObj, Method = "F-test", names = rnames, study = "Smoker")
lm <- readRDS("../analysis/data/raw_data/smoker_mlrt.rds")
lm <- augment(lm@qvalueObj, Method = "Moderated\nF-test", names = rnames, study = "Smoker")
l2 <- readRDS("../analysis/data/raw_data/smoker_lrt_boot.rds")
l2 <- augment(l2@qvalueObj, Method = "boot F-test", names = rnames, study = "Smoker")
lm2 <- readRDS("../analysis/data/raw_data/smoker_mlrt_boot.rds")
lm2 <- augment(lm2@qvalueObj, Method = "boot Moderated\nF-test", names = rnames, study = "Smoker")

df_smoker <- rbind(o, l, lm, l2, lm2)
saveRDS(df_smoker, file = "../analysis/data/smoker_comp.rds")


##### Dose study ####
gds <- getGEO("GSE4668") 
eset <- gds$GSE4668_series_matrix.txt.gz
rnames <- rownames(eset@assayData$exprs)
o <- readRDS("../analysis/data/raw_data/dose_odp.rds")
o <- augment(o@qvalueObj, Method = "ODP", names = rnames, study = "Dose")
l <- readRDS("../analysis/data/raw_data/dose_lrt.rds")
l <- augment(l@qvalueObj, Method = "F-test", names = rnames, study = "Dose")
lm <- readRDS("../analysis/data/raw_data/dose_mlrt.rds")
lm <- augment(lm@qvalueObj, Method = "Moderated\nF-test", names = rnames, study = "Dose")
l2 <- readRDS("../analysis/data/raw_data/dose_lrt_boot.rds")
l2 <- augment(l2@qvalueObj, Method = "boot F-test", names = rnames, study = "Dose")
lm2 <- readRDS("../analysis/data/raw_data/dose_mlrt_boot.rds")
lm2 <- augment(lm2@qvalueObj, Method = "boot Moderated\nF-test", names = rnames, study = "Dose")

df_dose <- rbind(o, l, lm, l2, lm2)
saveRDS(df_dose, file = "../analysis/data/dose_comp.rds")


#### Endotoxin ####
infile <- "../analysis/data/raw_data/"
rnames <- read.table(paste0(infile, "CR001_genes.txt"))$V1
o <- readRDS("../analysis/data/raw_data/endotoxin_odp.rds")
o <- augment(o@qvalueObj, Method = "ODP", names = rnames, study = "Endotoxin")
l <- readRDS("../analysis/data/raw_data/endotoxin_lrt.rds")
l <- augment(l@qvalueObj, Method = "F-test", names = rnames, study = "Endotoxin")
lm <- readRDS("../analysis/data/raw_data/endotoxin_mlrt.rds")
lm <- augment(lm@qvalueObj, Method = "Moderated\nF-test", names = rnames, study = "Endotoxin")
l2 <- readRDS("../analysis/data/raw_data/endotoxin_lrt_boot.rds")
l2 <- augment(l2@qvalueObj, Method = "boot F-test", names = rnames, study = "Endotoxin")
lm2 <- readRDS("../analysis/data/raw_data/endotoxin_mlrt_boot.rds")
lm2 <- augment(lm2@qvalueObj, Method = "boot Moderated\nF-test", names = rnames, study = "Endotoxin")

df_endotoxin <- rbind(o, l, lm, l2, lm2)
saveRDS(df_endotoxin, file = "../analysis/data/endotoxin_comp.rds")
