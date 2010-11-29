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
    pvals <- GSRI:::getPvalues(data, 1:nSamples, phenotype, test, testArgs)
    res <- GSRI:::getCDF(pvals, useGrenander=FALSE)
    if(useGrenander == TRUE)  {
      res$cdf <- GSRI:::cdfCorrect(res$sortedPval, res$cdf, 1-p)
      res$cdf <- GSRI:::grenanderInterp(res$sortedPvals, res$cdf)
    }
    nGenes <- length(pvals)
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
