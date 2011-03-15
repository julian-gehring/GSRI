rowt <- function(exprs, groups, id, testArgs)
  genefilter::rowttests(exprs[id, ], groups)$p.value


rowF <- function(exprs, groups, id, testArgs)
  genefilter::rowFtests(exprs[id, ], groups)$p.value


limmat <- function(exprs, groups, id, testArgs) {
  design <- cbind(offset=1, diff=groups)
  fit <- limma::lmFit(exprs, design)
  fit <- limma::eBayes(fit)
  pval <- fit$p.value[id,"diff"]

  return(pval)
}


limmatGsri <- function(exprs, groups, id, testArgs) {
  design <- cbind(offset=1, diff=groups)
  fit <- limma::lmFit(exprs, design)
  fit <- limma::eBayes(fit)
  pval <- fit$p.value[ ,"diff"]

  p0 <- les:::gsri(pval, grenander=testArgs$grenander, se=FALSE)[1]
  fit <- limma::eBayes(fit, p0)
  pval <- fit$p.value[id,"diff"]

  return(pval)
}
