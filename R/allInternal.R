calcGsri <- function(exprs, groups, name,
                     weight, grenander=TRUE, nBoot=100,
                     test="ttest", testArgs=NULL, alpha=0.05, ...) {

  nGenes <- nrow(exprs)
  if(ncol(exprs) != length(groups))
    stop("Number of columns of 'exprs' must match 'groups'.")
  if(!(is.function(test) && length(formals(test)) > 1))
    stop("'test' must be a function with at least two input arguments.")
  if(is.null(weight))
    weight <- rep(1, nGenes)

  pval <- multiStat(exprs, groups, weight, grenander, test, ...)
  les <- les:::fitGsri(pval, NULL, weight, nGenes, grenander, se=TRUE, custom=FALSE)
  p0 <- les[1]
  psd0 <- les[2]
  
  ## pval <- multiStat(exprs, 1:ncol(exprs), groups, test, testArgs)
  b <- boot::boot(t(exprs), gsriBoot, nBoot,
                  groups=groups, cweight=weight, grenander=grenander,
                  test=test, testArgs=testArgs)
  p1 <- max(b$t0[[1]], 0)
  bias <- apply(b$t, 2, mean) - b$t0
  pt <- b$t - bias
  gsri <- max(stats::quantile(pt, alpha, na.rm=TRUE), 0)
  
  p <- p1
  #pb <- replicate(nBoot,
  #                gsriBoot(exprs, groups, weight, grenander, test, ...))
  #p <- max(b$t0, 0)
  #psd <- stats::sd(pb)
  #gsri <- max(stats::quantile(pb, alpha, na.rm=TRUE), 0)
  nRegGenes <- as.integer(floor(p*nGenes))
  result <- data.frame(pRegGenes=p0, pRegGenesSd=psd0, nRegGenes=nRegGenes,
                       gsri=gsri, nGenes=nGenes, row.names=name)
  names(result)[4] <- sprintf("%s(%g%%)", "GSRI", alpha*100)
  res <- list(result=result, pval=pval)
  #res <- result
  
  return(res)
}


multiStat <- function(exprs, label, weight, grenander, test, ...) {

  pval <- test(exprs, label, ...)
  pval <- pval[!is.na(pval)] ## needed?

  return(pval)
}


gsriBoot <- function(exprs, index, groups, cweight, grenander, test, testArgs, ...) {

  pval <- multiStat(t(exprs), groups[index], cweight, grenander, test, ...)
  p <- les:::fitGsri(pval, NULL, cweight, length(pval), grenander, FALSE, FALSE)[1]

  return(p)
}


#gsriBoot <- function(exprs, groups, weight, grenander, test, ...) {
#
#  boot <- bootWithinGroup(exprs, groups)
#  res <- multiStat(boot$exprs, boot$label, weight, grenander, test, ...) ## weight[ord] !!!
#
#  return(res)
#}


bootWithinGroup <- function(exprs, label) {

  ord <- order(label)
  exprs <- exprs[ ,ord, drop=FALSE]
  label <- label[ord]

  nLevels <- nlevels(label) ## or length(nSamples) ? factor w/o group?
  nSamples <- tabulate(label)

  perm <- unlist(lapply(nSamples, sample.int, replace=TRUE))
  offset <- rep.int(cumsum(c(0, nSamples[-length(nSamples)])), nSamples)
  ordBoot <- perm + offset

  exprs <- exprs[ ,ordBoot, drop=FALSE]
  labelOrd <- label[ordBoot]
  if(!identical(label, labelOrd))
    warning("Bootstrapping lead to changes in the groups.")

  return(list(exprs=exprs, label=labelOrd))
}


getPvalues <-
function (exprs, d, groups, test, testArgs) 
{
    pvals <- GSRI:::multiStat(as.matrix(exprs), groups[d], test, 
        testArgs)
    pvals <- pvals[!is.na(pvals)]
    return(pvals)
}

