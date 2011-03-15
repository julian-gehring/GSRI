rowt <- function(exprs, groups, id, index, testArgs) {
  genefilter::rowttests(exprs[id,index,drop=FALSE], groups)$p.value
}


rowF <- function(exprs, groups, id, index, testArgs=list(var.equal=TRUE)) {
  genefilter::rowFtests(exprs[id,index,drop=FALSE], groups)$p.value
}


limmat <- function(exprs, groups, id, index,
                   testArgs=list(design=cbind(offset=1, diff=groups))) {
  fit <- limma::lmFit(exprs[ ,index,drop=FALSE], testArgs$design)
  fit <- limma::eBayes(fit)
  pval <- fit$p.value[id,"diff"]

  return(pval)
}

