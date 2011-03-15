getCDF <-
function (pvals, useGrenander) 
{
    nGenes <- length(pvals)
    ord <- sort.list(pvals, method = "quick", na.last = NA)
    pvals <- pvals[ord]
    uniquePvals <- unique(pvals)
    rcdf <- stats::ecdf(pvals)
    x <- environment(rcdf)$x
    cdf <- environment(rcdf)$y
    if (useGrenander == TRUE) 
        cdf <- GSRI:::grenanderInterp(x, cdf)
    if (nGenes != length(uniquePvals)) {
        tp <- as.numeric(table(pvals))
        w <- diff(c(0, cdf))/tp
        cdf <- rep.int(cdf, tp)
    }
    cdf <- cdf - 0.5/nGenes
    res <- list(sortedPvals = pvals, cdf = cdf)
    return(res)
}
getPvalues <-
function (data, d, phenotype, test, testArgs) 
{
    phenotype <- phenotype[d]
    pvals <- GSRI:::multiStat(as.matrix(data), phenotype, test, 
        testArgs)
    pvals <- pvals[!is.na(pvals)]
    return(pvals)
}
grenanderInterp <-
function (x, y) 
{
    ll <- gcmlcm(x, y, type = "lcm")
    pn <- vector("numeric", length(x))
    ind <- x %in% ll$x.knots
    pn[ind] <- ll$y.knots
    for (i in 2:length(ll$x.knots)) {
        xx <- ll$x.knots[(i - 1):i]
        yy <- ll$y.knots[(i - 1):i]
        ind <- (x > xx[1]) & (x < xx[2])
        xi <- x[ind]
        yi <- ll$slope.knots[i - 1] * (xi - xx[1]) + yy[1]
        pn[ind] <- yi
    }
    return(pn)
}
loadCls <-
function (fileName) 
{
    clsCont <- readLines(fileName)
    nSamples <- as.integer(unlist(strsplit(clsCont[[1]], " "))[1])
    phen <- unlist(strsplit(clsCont[[3]], " "))
    phenotype <- factor(phen)
    return(phenotype)
}
loadGct <-
function (fileName) 
{
    temp <- readLines(fileName, n = 3)
    colNames <- noquote(unlist(strsplit(temp[3], "\t")))[c(-1, 
        -2)]
    numCols <- length(colNames)
    colClasses <- c("character", "character", rep("numeric", 
        numCols))
    data <- utils::read.table(fileName, header = TRUE, skip = 2, row.names = 1, 
        na.strings = c("na", ""), sep = "\t", check.names = TRUE, 
        colClasses = colClasses, quote = "", as.is = TRUE)
    m <- as.matrix(data[, -1])
    return(m)
}
multiStat <-
function (data, phenotype, test, testArgs) 
{
    if (!is.function(test)) {
        if (test == "ttest") 
            pvals <- genefilter::rowttests(data, phenotype)$p.value
        if (test == "ftest") 
            pvals <- genefilter::rowFtests(data, phenotype)$p.value
    }
    else {
        if (is.null(testArgs)) 
            pvals <- test(data, phenotype)
        else pvals <- test(data, phenotype, testArgs)
    }
    return(pvals)
}
plotResults <-
function (res, geneSetName, p, gsri) 
{
    xfit <- c(0, 1)
    yfit <- p + (1 - p) * xfit
    graphics::plot(c(0, 1), c(p, p), col = "red", type = "l", 
        lty = 2, xlab = "p-values", ylab = "CDF(p)", main = geneSetName, 
        xlim = c(0, 1), ylim = c(0, 1))
    graphics::lines(c(0, 1), c(gsri[1], gsri[1]), col = "blue", 
        lty = 2)
    graphics::lines(xfit, yfit, col = "gray")
    graphics::points(res$sortedPvals, res$cdf, xlab = "p-values", 
        ylab = "CDF(p)", main = geneSetName, xlim = c(0, 1), 
        ylim = c(0, 1), pch = 20)
    graphics::text(1, p + 0.01, sprintf("%s=%.2f", "%RegGene", 
        p), cex = 0.8, adj = c(1, 0))
    graphics::text(1, gsri[1] - 0.01, sprintf("%s=%.2f", "GSRI", 
        gsri[1]), cex = 0.8, adj = c(1, 1))
}
slopeFast <-
function (x, y) 
{
    tx <- t(x)
    b <- as.numeric((tx %*% y)/(tx %*% x))
    return(b)
}
writeResults <-
function (res, geneSetName, p, prec) 
{
    xfit <- res$sortedPvals
    yfit <- p + (1 - p) * xfit
    dataFileName <- paste("GeneSet_", geneSetName, "_data.txt", 
        sep = "")
    dataResults <- list(sorted_pvals = signif(res$sortedPvals, 
        prec), cdf = signif(res$cdf, prec), fit = signif(yfit, 
        prec))
    utils::write.table(dataResults, file = dataFileName, quote = FALSE, 
        row.names = FALSE, sep = "\t")
}
gsriBoot <-
function (data, d, phenotype, useGrenander, test, testArgs) 
{
  pvals <- GSRI:::getPvalues(t(data), d, phenotype, test, testArgs)
  cdf <- GSRI:::getCDF(pvals, useGrenander=FALSE)
  res <- GSRI:::fitSlope(cdf$sortedPvals-1, cdf$cdf-1)
  if(useGrenander == TRUE)  {
    cdf$cdf <- GSRI:::cdfCorrect(cdf$sortedPval, cdf$cdf, 1-res)
    cdf$cdf <- GSRI:::grenanderInterp(cdf$sortedPvals, cdf$cdf)
    res <- GSRI:::fitSlope(cdf$sortedPvals-1, cdf$cdf-1)
  }
  return(res)
}
fitSlope <-
function(x, y)
{
  nValidGenes <- length(x)
  maxIter <- nValidGenes
  q <- 1
  restOld <- 0
  for (nIterate in 1:maxIter) {
    rest <- nValidGenes - ceiling(q * nValidGenes)
    rest <- max(c(restOld, rest, 1))
    rest <- min(nValidGenes - 1, rest)
    if (is.na(rest) || restOld == rest) 
      break
    ind <- rest:nValidGenes
    q <- GSRI:::slopeFast(x[ind], y[ind])
    restOld <- rest
  }
  res <- 1 - q
  return(res)
}
cdfCorrect <-
function(x, y, q0)
{
  z <- y
  indLower <- y < q0*x
  indUpper <- y > 1 - q0*(1 - x)
  z[indLower] <- q0*x[indLower]
  z[indUpper] <- 1 - q0*(1 - x[indUpper])
  return(z)
}
calcGsri <- function(...) {
  ## dummy
}
