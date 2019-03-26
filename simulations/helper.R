###############################################################################
# Calculate the row variances
# x: matrix
###############################################################################
rowVars <- function(x) {
  rowSums((x - rowMeans(x)) ^ 2) / (dim(x)[2] - 1)
}

###############################################################################
# Calculate the Signal-to-Noise Ratio for each gene
# x: deObj from edge
# fdr_cutoff: false discovery rate cutoff
###############################################################################
get_SNR <- function(x, weights = NULL, fdr_cutoff = 0.05) {
  sig_id <- which(edge::qvalueObj(x)$qvalues < fdr_cutoff)
  null_id <- which(edge::qvalueObj(x)$qvalues > fdr_cutoff)
  if (!is.null(weights)) {
    edgeFit <- edge::fit_models(x, weights = weights)
  } else {
    edgeFit <- edge::fit_models(x)
  }
  sigFit <- edgeFit@fit.full[sig_id,]
  resFit <- edgeFit@res.full[sig_id,]
  nullFit <- edgeFit@fit.null[null_id,]
  nresFit <- edgeFit@res.null[null_id,]
  var_signal <- c(rowVars(sigFit), rowVars(nullFit))
  var_noise <- c(rowVars(resFit), rowVars(nresFit))
  return(list(snr = var_signal / var_noise, alt_id = sig_id, null_id = null_id))
}

###############################################################################
# Calculate estimated FDR
# oracle: boolean vector indicating which genes are DE
# pvals: vector of p-values
# fdr: vector of false discovery rate cutoffs
###############################################################################
estFDR <- function(oracle, pvals, fdr = seq(0, .15, 0.001)) {
  qvals <- qvalue::qvalue(pvals)$qvalue
  est_fdr <- tot_sigGenes <- rep(0, length(fdr))
  for (i in 1:length(fdr)) {
    est <- oracle[qvals < fdr[i]]
    tot_sigGenes[i] <- length(est)
    FP = sum(!est)
    est_fdr[i] <- FP / pmax(tot_sigGenes[i], 1)
  }
  return(data.frame(true_fdr = fdr,
                    tot_genes = tot_sigGenes,
                    est_fdr = est_fdr))
}
