calcGsri <- function(exprs, groups, name, id, weights,
                     grenander=TRUE, nBoot=100,
                     test=NULL, testArgs=NULL, alpha=0.05, ...) {
  
  if(length(id) == 0 || !is.logical(id))
    stop("Gene set has no matches in the expression data.")
  if(ncol(exprs) != length(groups))
    stop("Number of columns of 'exprs' must match 'groups'.")
  if(!(is.function(test) && length(formals(test)) > 1))
    stop("'test' must be a function with at least two input arguments.")
  if(length(alpha) > 1) {
    alpha <- alpha[1]
    warning("Taking only first value of 'alpha'.")
  }
  nGenesGs <- sum(id)
  if(nGenesGs == 0) {
    result <- data.frame(pRegGenes=NA, pRegGenesSd=NA, nRegGenes=NA,
                         gsri=NA, nGenes=nGenesGs, row.names=name)
    res <- list(result=result, pval=NULL)
    return(res)
  }

  ## reorder data for bootstrapping
  ord <- order(groups)
  exprs <- exprs[ ,ord]
  groups <- groups[ord]
  nSamples <- tabulate(groups)

  ## calculate Lambda
  pval <- multiStat(exprs, groups, id, 1:ncol(exprs), test, testArgs)
  nPval <- length(pval)
  nGenes <- nrow(exprs)
  if(is.null(weights))
    weights <- rep(1/nGenes, nGenes)
  if(length(weights) == nGenes)
    weights <- weights[id]
  if(length(weights) != nGenesGs)
    stop("'weight' must contain one value for each gene.")
  if(length(pval) != nGenesGs)
    stop("'test' must return one p-value for each gene in the gene set.")
  cdf <- les:::wcdfGrenander(pval, weights, nPval, grenander, FALSE)
  l0 <- les:::itLinReg(cdf$pval, cdf$cdf, weights, nPval, FALSE, FALSE, FALSE)

  ## bootstrapping
  lb <- replicate(nBoot,
                  gsriBoot(exprs, groups, weights, id, grenander, test, testArgs, nSamples))
  lb <- lb - mean(lb) + l0
  lsd <- stats::sd(lb)
  gsri <- max(stats::quantile(lb, alpha, na.rm=TRUE), 0)
  nRegGenes <- as.integer(floor(l0*nGenesGs))

  ## gsri results
  result <- data.frame(pRegGenes=l0, pRegGenesSd=lsd, nRegGenes=nRegGenes,
                       gsri=gsri, nGenes=nGenesGs, row.names=name)
  names(result)[4] <- sprintf("%s(%g%%)", "GSRI", alpha*100)

  ## cdf results
  geneNames <- rownames(exprs[id, ])[order(pval)]
  wcdf <- data.frame(pval=cdf$pval, cdf=cdf$cdf, row.names=geneNames)
  
  res <- list(result=result, cdf=wcdf)

  return(res)
}


multiStat <- function(exprs, groups, id, index, test, testArgs=NULL) {

  if(is.null(testArgs))
    pval <- test(exprs, groups, id, index)
  else
    pval <- test(exprs, groups, id, index, testArgs)
  pval <- pval[!is.na(pval)]
  if(length(pval) != sum(id))
    stop("Test statistics must return one p-value for each gene in the gene set.")

  return(pval)
}


gsriBoot <- function(exprs, groups, weights, id, grenander, test, testArgs, nSamples) {

  index <- bootInGroups(nSamples)
  pval <- multiStat(exprs, groups, id, index, test, testArgs)
  res <- les:::fitGsri(pval, NULL, weights, length(pval), grenander, se=FALSE, custom=FALSE)[1]

  return(res)
}


bootInGroups <- function(nSamples) {

  perm <- unlist(lapply(nSamples, sample.int, replace=TRUE))
  offset <- rep.int(cumsum(c(0L, nSamples[-length(nSamples)])), nSamples)
  index <- perm + offset
  ## check: only one group, groups with one sample, etc.

  return(index)
}


getArgs <- function(name, first=NULL, last=NULL, ...) {

  ind <- which(names(...) %in% name)[1] ## [[]] allows only one element for indexing
  middle <- if(length(ind) != 0) ...[[ind]] else NULL
  args <- c(first, middle, last)
  args <- args[!duplicated(names(args))]
  
  return(args)
}

