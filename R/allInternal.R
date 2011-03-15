calcGsri <- function(exprs, groups, name, id, 
                     weight, grenander=TRUE, nBoot=100,
                     test="ttest", testArgs=NULL, alpha=0.05, ...) {

  if(length(id) == 0 || !all(is.logical(id)))
    stop("Gene set has no matches in the expression data.")
  if(ncol(exprs) != length(groups))
    stop("Number of columns of 'exprs' must match 'groups'.")
  if(!(is.function(test) && length(formals(test)) > 1))
    stop("'test' must be a function with at least two input arguments.")
  nGenesGs <- sum(id)
  
  pval <- multiStat(exprs, groups, id, grenander, test, testArgs)
  nPval <- length(pval)
  if(is.null(weight))
    weight <- rep(1, nPval)
  les <- les:::fitGsri(pval, NULL, weight, nPval, grenander, se=TRUE, custom=FALSE)
  p0 <- les[1]
  psd0 <- les[2]
  
  ## pval <- multiStat(exprs, 1:ncol(exprs), groups, test, testArgs)
  b <- boot::boot(t(exprs), gsriBoot, nBoot,
                  groups=groups, id=id, cweight=weight, 
                  grenander=grenander, test=test, testArgs=testArgs)
  p1 <- max(b$t0[[1]], 0)
  bias <- apply(b$t, 2, mean) - b$t0
  pt <- b$t - bias
  gsri <- max(stats::quantile(pt, alpha, na.rm=TRUE), 0)
  
  p <- p1
  #pb <- replicate(nBoot,
  #                gsriBoot(exprs, groups, weight, grenander, test, testArgs))
  #p <- max(b$t0, 0)
  #psd <- stats::sd(pb)
  #gsri <- max(stats::quantile(pb, alpha, na.rm=TRUE), 0)
  nRegGenes <- as.integer(floor(p*nGenesGs))
  result <- data.frame(pRegGenes=p0, pRegGenesSd=psd0, nRegGenes=nRegGenes,
                       gsri=gsri, nGenes=nGenesGs, row.names=name)
  names(result)[4] <- sprintf("%s(%g%%)", "GSRI", alpha*100)
  res <- list(result=result, pval=pval)
  #res <- result
  
  return(res)
}


multiStat <- function(exprs, groups, id, grenander, test, testArgs) {

  pval <- test(exprs, groups, id, testArgs)
  pval <- pval[!is.na(pval)] ## needed?

  return(pval)
}


gsriBoot <- function(exprs, index, groups, id, cweight, grenander, test, testArgs, ...) {

  pval <- multiStat(t(exprs), groups[index], id, grenander, test, testArgs)
  p <- les:::fitGsri(pval, NULL, cweight, length(pval), grenander, FALSE, FALSE)[1]

  return(p)
}


#gsriBoot <- function(exprs, groups, weight, grenander, test, testArgs) {
#
#  boot <- bootWithinGroup(exprs, groups)
#  res <- multiStat(boot$exprs, boot$label, weight, grenander, test, testArgs) ## weight[ord] !!!
#
#  return(res)
#}


bootWithinGroup <- function(exprs, groups) {

  ord <- order(groups)
  exprs <- exprs[ ,ord, drop=FALSE]
  groups <- groups[ord]

  nLevels <- nlevels(groups) ## or length(nSamples) ? factor w/o group?
  nSamples <- tabulate(groups)

  perm <- unlist(lapply(nSamples, sample.int, replace=TRUE))
  offset <- rep.int(cumsum(c(0, nSamples[-length(nSamples)])), nSamples)
  ordBoot <- perm + offset

  exprs <- exprs[ ,ordBoot, drop=FALSE]
  groupsOrd <- groups[ordBoot]
  if(!identical(groups, groupsOrd))
    warning("Bootstrapping lead to changes in the groups.")

  return(list(exprs=exprs, groups=groupsOrd))
}


