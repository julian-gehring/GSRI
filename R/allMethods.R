## gsri ##
setGeneric("gsri",
           function(exprs, groups, geneSet, names=NULL, weight=NULL, nBoot=100, 
                   test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE, ...)
           standardGeneric("gsri"))

setMethod("gsri",
          signature("matrix", "factor", "missing"),
          function(exprs, groups, geneSet, names=NULL, weight=NULL, nBoot=100, 
                   test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE,
                   id=!logical(nrow(exprs)), ...) {
            
            res <- calcGsri(exprs, groups, names, id, weight,
                            grenander, nBoot, test, testArgs, alpha)
            cdf <- list(res$cdf)
            names(cdf) <- names
            parms <- list(weight=weight, nBoot=nBoot, test=test, alpha=alpha,
                          grenander=grenander, testArgs=testArgs)
            
            
            object <- new("Gsri",
                          result=res$result, cdf=cdf, parms=parms)
            
            return(object)
          })

setMethod("gsri",
          signature("ExpressionSet", "factor", "missing"),
          function(exprs, groups, geneSet, names=NULL, weight=NULL, nBoot=100, 
                   test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE, ...) {

            object <- gsri(exprs(exprs), groups, names=names, weight=weight,
                           nBoot=nBoot, test=test, testArgs=testArgs, alpha=alpha,
                           grenander=grenander)
            
            return(object)
          })

setMethod("gsri",
          signature("matrix", "factor", "GeneSet"),
          function(exprs, groups, geneSet, names=NULL, weight=NULL, nBoot=100, 
                   test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE, ...) {

            if(is.null(names))
              names <- setName(geneSet)
            if(is.na(names))
              names <- NULL
            id <-  rownames(exprs) %in% geneIds(geneSet)
            object <- gsri(exprs, groups, names=names, weight=weight,
                           nBoot=nBoot, test=test, testArgs=testArgs, alpha=alpha,
                           grenander=grenander, id=id)

            return(object)
          })

setMethod("gsri",
          signature("ExpressionSet", "factor", "GeneSet"),
          function(exprs, groups, geneSet, names=NULL, weight=NULL, nBoot=100, 
                   test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE, ...) {

            object <- gsri(exprs(exprs), groups, geneSet, names=names, weight=weight,
                           nBoot=nBoot, test=test, testArgs=testArgs, alpha=alpha,
                           grenander=grenander)
            
            return(object)
          })

setMethod("gsri",
          signature("matrix", "factor", "GeneSetCollection"),
          function(exprs, groups, geneSet, names=NULL, weight=NULL, nBoot=100, 
                   test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE, nCores=NULL, ...) {

            if(is.null(names))
              names <- names(geneSet)
            res <- les:::mcsapply(geneSet, gsri, exprs=exprs, groups=groups, name=NULL,
                                  weight=weight, nBoot=nBoot, grenander=grenander, test=test,
                                  testArgs=testArgs, alpha=alpha, mc.cores=nCores)

            object <- new("Gsri",
                          result=as.data.frame(do.call(rbind, lapply(res, getGsri)), row.names=names),
                          cdf=sapply(res, getCdf),
                          parms=getParms(res[[1]])
                          )

            return(object)
          })

setMethod("gsri",
          signature("ExpressionSet", "factor", "GeneSetCollection"),
          function(exprs, groups, geneSet, names=NULL, weight=NULL, nBoot=100, 
                   test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE, nCores=NULL, ...) {

            object <- gsri(exprs(exprs), groups, geneSet, names=names, weight=weight,
                         nBoot=nBoot, test=test, testArgs=testArgs, alpha=alpha,
                         grenander=grenander, nCores=nCores)
            
            return(object)
          })


## getGsri ##
setGeneric("getGsri",
           function(object, index, ...)
           standardGeneric("getGsri"))

setMethod("getGsri",
          signature("Gsri", "missing"),
          function(object, index) {
            return(object@result)
          })

setMethod("getGsri",
          signature("Gsri", "ANY"),
          function(object, index) {
            return(object@result[index, ])
          })


## getCdf ##
setGeneric("getCdf",
           function(object, index, ...)
           standardGeneric("getCdf"))

setMethod("getCdf",
          signature("Gsri", "missing"),
          function(object, index) {
            res <- object@cdf
#            if(length(res) == 1)
#              res <- res[[1]]
            return(res)
          })

setMethod("getCdf",
          signature("Gsri", "ANY"),
          function(object, index) {
            return(object@cdf[[index]])
          })


## getParms ##
setGeneric("getParms",
           function(object, ...)
           standardGeneric("getParms"))

setMethod("getParms",
          signature("Gsri"),
          function(object) {
            return(object@parms)
          })


## export ##
setGeneric("export",
           function(object, file, ...)
           standardGeneric("export"))

setMethod("export",
          signature("Gsri", "character"),
          function(object, file, digits=Inf) {
            result <- round(getGsri(object), digits)
            write.table(result, file, sep="\t")
          })


## show ##
setMethod("show",
          signature("Gsri"),
          function(object) {
            print(object@result)
          })


## summary ##
setGeneric("summary",
           function(object, ...)
           standardGeneric("summary"))

setMethod("summary",
          signature("Gsri"),
          function(object, ...) {
            show(object)
          })


## sortGsri ##
setGeneric("sortGsri",
           function(x, names, decreasing=TRUE, na.last=NA, ...)
           standardGeneric("sortGsri"))

setMethod("sortGsri",
          signature("Gsri"),
          function(x, names, decreasing=TRUE, na.last=NA) {

            res <- getGsri(x)
            sub <- subset(res, select=names)
            ord <- do.call(order, c(sub, list(decreasing=decreasing, na.last=na.last)))
            x@result <- res[ord, ]
            x@cdf <- x@cdf[ord]

            return(x)
          })


## plot ##
setGeneric("plot",
           function(x, y, ...)
           standardGeneric("plot"))

setMethod("plot",
          signature("Gsri", "ANY"),
          function(x, y, digits=2, ...) {
            
            result <- getGsri(x)
            sel <- result[y, ]
            if(nrow(sel) > 1) {
              sel <- sel[1, ]
              warning("More than one gene set chosen, taking only the first one.")
            }
            ind <- which(rownames(result) %in% rownames(sel))
            if(all(is.na(sel)) || length(ind) == 0)
              stop("No valid index for selecting a gene set.")
            plot(x, ind, digits, ...)
          })

setMethod("plot",
          signature("Gsri", "integer"),
          function(x, y, digits=2, ...) {

            args <- list(...)
            result <- getGsri(x)[y, ]
            
            p <- as.numeric(result[ ,1])  ## do with names
            g <- as.numeric(result[ ,4])  ## do with names

            ## plot arguments
            ## plot
            plot1 <- list(x=NA, y=NA)
            plot3 <- list(xlab="p-values", ylab="ECDF(p)", main=rownames(result),
                          xlim=c(0, 1), ylim=c(0, 1))
            plotArgs <- getArgs("plot", plot1, plot3, args)
            ## reg
            reg1 <- list(x=c(0, 1), y=c(p, p))
            reg3 <- list(col="red", type="l", lty=2)
            regArgs <- getArgs("reg", reg1, reg3, args)
            ## gsri
            gsri1 <- list(x=c(0, 1), y=c(g, g))
            gsri3 <- list(col="blue", lty=2)
            gsriArgs <- getArgs("gsri", gsri1, gsri3, args)
            ## fit
            fit1 <- list(x=c(0, 1), y=p+(1-p)*c(0, 1))
            fit3 <- list(col="gray")
            fitArgs <- getArgs("fit", fit1, fit3, args)
            ## ecdf
            ecdf1 <- list(x=x@cdf[[y]]$pval, y=x@cdf[[y]]$cdf)
            ecdf3 <- list(type="p", pch=20)
            ecdfArgs <- getArgs("ecdf", ecdf1, ecdf3, args)

            format <- paste("%s=%.", digits, "f", sep="")
            ## regText
            regText1 <- list(x=1, y=p+0.01)
            regText3 <- list(labels=sprintf(format, "%RegGene", p), cex=0.8, adj=c(1, 0))
            regTextArgs <- getArgs("regText", regText1, regText3, args)
            ## gsriText
            gsriText1 <- list(x=1, y=g-0.01)
            gsriText3 <- list(labels=sprintf(format, "%GSRI", g), cex=0.8, adj=c(1, 1))
            gsriTextArgs <- getArgs("gsriText", gsriText1, gsriText3, args)
            
            ## plot calls
            do.call("plot", plotArgs)
            do.call("lines", regArgs)
            do.call("lines", gsriArgs)
            do.call("lines", fitArgs)
            do.call("lines", ecdfArgs)
            do.call("text", regTextArgs)
            do.call("text", gsriTextArgs)
          })


## readCls ##
setGeneric("readCls",
           function(file, ...)
           standardGeneric("readCls"))

setMethod("readCls",
          signature("character"),
          function(file, ...) {
            if(!file.exists(file))
              stop(sprintf("%s '%s' %s", "File", file, "does not exist."))
            clsCont <- readLines(file)
            header <- as.integer(unlist(strsplit(clsCont[[1]], " ")))
            groups <- factor(unlist(strsplit(clsCont[[3]], " ")))
            if(length(groups) != header[1] || nlevels(groups) != header[2])
              warning(sprintf("%s '%s' %s", "Data in file", basename(file),
                              "is not consistent."))
            
            return(groups)
          })


## readGct ##
setGeneric("readGct",
           function(file, ...)
           standardGeneric("readGct"))

setMethod("readGct",
          signature("character"),
          function(file, ...) {
            if(!file.exists(file))
              stop(sprintf("%s '%s' %s", "File", file, "does not exist."))
            header <- readLines(file, n=3)
            extend <- as.integer(noquote(unlist(strsplit(header[2], "\t")))[c(1,2)])
            exprs <- utils::read.table(file, header=TRUE, skip=2, as.is=TRUE,
                                      row.names=1, sep="\t", quote="",
                                      na.strings=c("na", ""))
            exprs <- as.matrix(exprs[ ,-1])
            if(any(dim(exprs) != extend))
              warning(sprintf("%s '%s' %s", "Data in file", basename(file),
                              "is not consistent."))
            
            return(exprs)
          })
