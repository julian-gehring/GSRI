getCDF <-
function (pvals, useGrenander) 
{
    nGenes <- length(pvals)
    cor <- 0.5/nGenes
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
    cdf <- cdf - cor
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
    ll <- fdrtool::gcmlcm(x, y, type = "lcm")
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
gsri <-
function (data, phenotype, geneSetName, useGrenander = FALSE, 
    plotResults = TRUE, writeResults = FALSE, nBootstraps = 100, 
    test = "ttest", testArgs = NULL, alpha = 0.05, prec = 4) 
{
    nSamples <- ncol(data)
    if (!is.function(test) && !is.character(test)) 
        stop("Argument 'test' must be of class 'function' or 'character'")
    if (is.character(test) && is.na(pmatch(test, c("ttest", "ftest"), 
        NA))) 
        stop("Character string for 'test' must be 'ttest' or 'ftest'")
    if (nSamples != length(phenotype)) 
        stop("Sizes of 'data' and 'phenotype' do not match.")
    if (!is.factor(phenotype)) 
        phenotype <- factor(phenotype)
    b <- boot::boot(t(data), GSRI:::gsriBoot, nBootstraps, phenotype = phenotype, 
        useGrenander = useGrenander, test = test, testArgs = testArgs)
    p <- max(b$t0[[1]], 0)
    psd <- stats::sd(b$t[, 1])
    bias <- apply(b$t, 2, mean) - b$t0
    pt <- b$t - bias
    gsri <- max(stats::quantile(pt, alpha, na.rm = TRUE), 0)
    pvals <- GSRI:::getPvalues(data, 1:nSamples, phenotype, test, 
        testArgs)
    nGenes <- length(pvals)
    res <- GSRI:::getCDF(pvals, useGrenander)
    numRegGenes <- p * nGenes
    numRegGenesSd <- psd * nGenes
    if (plotResults == TRUE) 
        GSRI:::plotResults(res, geneSetName, p, gsri)
    if (writeResults == TRUE) 
        GSRI:::writeResults(res, geneSetName, p, prec)
    result <- list(geneSet = geneSetName, percRegGenes = p, percRegGenesSd = psd, 
        numRegGenes = numRegGenes, numRegGenesSd = numRegGenesSd, 
        gsri = gsri, nGenes = nGenes)
    return(result)
}
gsriBoot <-
function (inputData, d, phenotype, useGrenander, test, testArgs) 
{
    data <- t(inputData)
    pvals <- GSRI:::getPvalues(data, d, phenotype, test, testArgs)
    nValidGenes <- length(pvals)
    maxIter <- nValidGenes + 2
    p <- 1
    restOld <- 0
    res <- GSRI:::getCDF(pvals, useGrenander)
    for (nIterate in 1:maxIter) {
        rest <- nValidGenes - ceiling(p * nValidGenes)
        rest <- max(c(restOld, rest, 1))
        rest <- min(nValidGenes - 1, rest)
        if (restOld == rest) 
            break
        x <- res$sortedPvals[rest:nValidGenes] - 1
        y <- res$cdf[rest:nValidGenes] - 1
        p <- GSRI:::slopeFast(x, y)
        restOld <- rest
        if (is.na(p)) 
            break
    }
    erg <- 1 - p
    return(erg)
}
gsriFromFile <-
function (dataFileName, phenotypeFileName, geneSetFileName, useGrenander = FALSE, 
    plotResults = FALSE, writeResults = FALSE, writeSummary = FALSE, 
    minGeneSetSize = 10, nBootstraps = 100, test = "ttest", testArgs = NULL, 
    alpha = 0.05, verbose = TRUE, prec = 4) 
{
    if (file.exists(dataFileName) == FALSE) 
        stop(sprintf("%s '%s' %s", "File", dataFileName, "cannot be read"))
    if (file.exists(phenotypeFileName) == FALSE) 
        stop(sprintf("%s '%s' %s", "File", phenotypeFileName, 
            "cannot be read"))
    if (file.exists(geneSetFileName) == FALSE) 
        stop(sprintf("%s '%s' %s", "File", geneSetFileName, "cannot be read"))
    data <- GSRI:::loadGct(dataFileName)
    if (verbose == TRUE) 
        message(sprintf("%s '%s' %s", "Data file", basename(dataFileName), 
            "loaded"))
    phenotype <- GSRI:::loadCls(phenotypeFileName)
    if (verbose == TRUE) 
        message(sprintf("%s '%s' %s", "Phenotype file",
                        basename(phenotypeFileName), "loaded"))
    geneSets <- readLines(geneSetFileName)
    res <- list()
    c <- 1
    for (i in 1:length(geneSets)) {
        temp <- strsplit(geneSets[i], "\t")[[1]]
        geneSetName <- temp[1]
        gsData <- data[intersect(row.names(data), temp[3:length(temp)]), 
            ]
        gsSetSize <- nrow(gsData)
        gsGenesResults <- paste(row.names(gsData), collapse = " : ")
        if (gsSetSize >= minGeneSetSize) {
            if (verbose == TRUE) 
                message(sprintf("%s '%s' (%d %s)", "Calculating GSRI for", 
                  geneSetName, gsSetSize, "genes"))
            gsriRes <- GSRI::gsri(gsData, phenotype, geneSetName, 
                useGrenander = useGrenander, plotResults = plotResults, 
                writeResults = writeResults, test = test, testArgs = testArgs, 
                alpha = alpha, nBootstraps = nBootstraps, prec = prec)
            res$geneSet[c] <- gsriRes$geneSet
            res$percRegGenes[c] <- signif(gsriRes$percRegGenes, 
                prec)
            res$percRegGenesSd[c] <- signif(gsriRes$percRegGenesSd, 
                prec)
            res$numRegGenes[c] <- signif(gsriRes$numRegGenes, 
                prec)
            res$numRegGenesSd[c] <- signif(gsriRes$numRegGenesSd, 
                prec)
            res$gsri[c] <- signif(gsriRes$gsri, prec)
            res$nGenes[c] <- signif(gsriRes$nGenes, prec)
            res$geneSetGenes[c] <- gsGenesResults
            c <- c + 1
        }
    }
    if (writeSummary == TRUE) {
        resultFileName <- paste(strsplit(dataFileName, ".gct")[[1]], 
            "_results.txt", sep = "")
        resultFile <- file.path(getwd(), resultFileName)
        write.table(res, file = resultFile, quote = TRUE, row.names = FALSE, 
            sep = "\t")
    }
    return(res)
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
