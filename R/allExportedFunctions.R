rowt <- function(exprs, groups, id, index, testArgs) {
  #genefilter:::rowcoltt(exprs[id, ,drop=FALSE], groups, FALSE, 1L)$p.value
  genefilter::rowttests(exprs[id,index,drop=FALSE], groups)$p.value
}


rowF <- function(exprs, groups, id, index, testArgs=TRUE) {
  #genefilter:::rowcolFt(exprs[id, ,drop=FALSE], groups, testArgs, 1L)$p.value
  genefilter::rowFtests(exprs[id,index,drop=FALSE], groups)$p.value
}


limmat <- function(exprs, groups, id, index, testArgs) {
  design <- cbind(offset=1, diff=groups)
  fit <- limma::lmFit(exprs[ ,index,drop=FALSE], design)
  fit <- limma::eBayes(fit)
  pval <- fit$p.value[id,"diff"]

  return(pval)
}


limmatGsri <- function(exprs, groups, id, index, testArgs) {
  design <- cbind(offset=1, diff=groups)
  fit <- limma::lmFit(exprs[ ,index,drop=FALSE], design)
  fit <- limma::eBayes(fit)
  pval <- fit$p.value[ ,"diff"]

  p0 <- les:::gsri(pval, grenander=testArgs$grenander, se=FALSE)[1]
  fit <- limma::eBayes(fit, p0)
  pval <- fit$p.value[id,"diff"]

  return(pval)
}
