library(splines)
library(tidyverse)

expression <- read_delim("../data/raw_data/kidney_expr.txt", delim = "\t", col_names = FALSE)
cov <-  t(read_delim("../data/raw_data/kidney_cov.txt", delim = "\t", col_names = FALSE))
colnames(cov) <- c("tissue", "age","pathology", "creatinine", "sex","gfr")
cov = type_convert(as_tibble(cov[-1,]))
expression <- log2(expression + 10)
rnames <- read_delim("../data/raw_data/kidney_genes.txt", delim = "\t", col_names = FALSE)

resid = expression[, 1:72]
ID <- rep(FALSE, nrow(expression))
for (i in 1:nrow(expression)) {
  print(i)
  for (j in 1:72){
    s <- cov$sex[-j]
    a = cov$age[-j]
    b = cov$tissue[-j]
    out <- lm( as.numeric(expression[i,-j])~s + b + ns(a, df = 9))
    resid[i,j] <- expression[i,j] - predict(out, data.frame(a = cov$age[j], b = cov$tissue[j], s =cov$sex[j]))
  }
}
saveRDS(resid, file = "../data/raw_data/kidney_residuals.rds")
filter.probes <- apply(resid, 1, FUN = function(x) (min(x) < -3 * IQR(x) + median(x)) | (max(x) > 3 * IQR(x) + median(x)))

saveRDS(!filter.probes, file = "../data/filter_probes.rds")
