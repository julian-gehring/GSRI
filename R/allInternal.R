calcGsri <- function(data, phenotype, name,
                     weight, grenander=TRUE, nBoot=100,
                     test="ttest", testArgs=NULL, alpha=0.05, ...) {

  nGenes <- nrow(data)
  if(ncol(data) != length(phenotype))
    stop("Number of columns of 'data' must match 'phenotype'.")
  if(!(is.function(test) && length(formals(test)) > 1))
    stop("'test' must be a function with at least two input arguments.")
  if(is.null(weight))
    weight <- rep(1, nGenes)

  pval <- multiStat(data, phenotype, weight, grenander, test, ...)
  les <- les:::fitGsri(pval, NULL, weight, nGenes, grenander, se=TRUE, custom=FALSE)
  p0 <- les[1]
  psd0 <- les[2]
  
  ## pval <- multiStat(data, 1:ncol(data), phenotype, test, testArgs)
  b <- boot::boot(t(data), gsriBoot, nBoot,
                  phenotype=phenotype, cweight=weight, grenander=grenander,
                  test=test, testArgs=testArgs)
  p1 <- max(b$t0[[1]], 0)
  bias <- apply(b$t, 2, mean) - b$t0
  pt <- b$t - bias
  gsri <- max(stats::quantile(pt, alpha, na.rm=TRUE), 0)
  
  p <- p1
  #pb <- replicate(nBoot,
  #                gsriBoot(data, phenotype, weight, grenander, test, ...))
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


multiStat <- function(data, label, weight, grenander, test, ...) {

  pval <- test(data, label, ...)
  pval <- pval[!is.na(pval)] ## needed?

  return(pval)
}


gsriBoot <- function(data, index, phenotype, cweight, grenander, test, testArgs, ...) {

  pval <- multiStat(t(data), phenotype[index], cweight, grenander, test, ...)
  p <- les:::fitGsri(pval, NULL, cweight, length(pval), grenander, FALSE, FALSE)[1]

  return(p)
}


#gsriBoot <- function(data, phenotype, weight, grenander, test, ...) {
#
#  boot <- bootWithinGroup(data, phenotype)
#  res <- multiStat(boot$data, boot$label, weight, grenander, test, ...) ## weight[ord] !!!
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
    warning("Bootstrapping lead to changes in the phenotype.")

  return(list(data=exprs, label=labelOrd))
}


getPvalues <-
function (data, d, phenotype, test, testArgs) 
{
    pvals <- GSRI:::multiStat(as.matrix(data), phenotype[d], test, 
        testArgs)
    pvals <- pvals[!is.na(pvals)]
    return(pvals)
}

