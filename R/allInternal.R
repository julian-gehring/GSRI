calcGsri <- function(exprs, groups, name, id, weights,
                     grenander=TRUE, nBoot=100,
                     test="ttest", testArgs=NULL, alpha=0.05, ...) {

  if(length(id) == 0 || !is.logical(id))
    stop("Gene set has no matches in the expression data.")
  if(ncol(exprs) != length(groups))
    stop("Number of columns of 'exprs' must match 'groups'.")
  if(!(is.function(test) && length(formals(test)) > 1))
    stop("'test' must be a function with at least two input arguments.")
  nGenesGs <- sum(id)
  
  pval <- multiStat(exprs, groups, id, grenander, test, testArgs)
  nPval <- length(pval)
  if(is.null(weights))
    weights <- rep(1/nGenesGs, nGenesGs)
  if((length(weights) != nGenesGs) || (length(weights) != nPval))
    stop("'weight' must contain one value for each gene in the gene set.")
  les <- les:::fitGsri(pval, NULL, weights, nPval, grenander, se=FALSE, custom=FALSE)
  p <- les[1]
  
  pb <- replicate(nBoot, gsriBoot(exprs, groups, weights, id, grenander, test, testArgs))
  pb <- pb - mean(pb) + p
  psd <- stats::sd(pb)
  gsri <- max(stats::quantile(pb, alpha, na.rm=TRUE), 0)

  nRegGenes <- as.integer(floor(p*nGenesGs))
  result <- data.frame(pRegGenes=p, pRegGenesSd=psd, nRegGenes=nRegGenes,
                       gsri=gsri, nGenes=nGenesGs, row.names=name)
  names(result)[4] <- sprintf("%s(%g%%)", "GSRI", alpha*100)
  res <- list(result=result, pval=pval)

  return(res)
}


multiStat <- function(exprs, groups, id, grenander, test, testArgs) {

  pval <- test(exprs, groups, id, testArgs)
  if(length(pval) != sum(id))
    stop("Test statistics must return one p-value for each gene in the gene set.")
  #pval <- pval[!is.na(pval)] ## needed?

  return(pval)
}


gsriBoot <- function(exprs, groups, weights, id, grenander, test, testArgs) {

  bo <- bootWithinGroup(exprs, groups)
  pval <- multiStat(bo$exprs, bo$groups, id, grenander, test, testArgs)
  res <- les:::fitGsri(pval, NULL, weights, length(pval), grenander, se=TRUE, custom=FALSE)[1]

  return(res)
}


bootWithinGroup <- function(exprs, groups, weights) {

  ord <- order(groups)
  exprs <- exprs[ ,ord,drop=FALSE]
  groups <- groups[ord] ## needed for correct resampling

  nLevels <- nlevels(groups) ## or length(nSamples) ? factor w/o group?
  nSamples <- tabulate(groups)

  perm <- unlist(lapply(nSamples, sample.int, replace=TRUE))
  offset <- rep.int(cumsum(c(0, nSamples[-length(nSamples)])), nSamples)
  ordBoot <- perm + offset

  exprs <- exprs[ ,ordBoot,drop=FALSE]
  if(!identical(groups, groups[ordBoot])) ## remove? just control
    warning("Bootstrapping lead to changes in the groups.")
  res <- list(exprs=exprs, groups=groups)

  return(res)
}

