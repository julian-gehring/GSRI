## gsri ##
setGeneric("gsri",
           function(exprs, phenotype, geneSet, names=NULL, weight=NULL, nBoot=100, 
                   test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE)
           standardGeneric("gsri"))

setMethod("gsri",
          signature("matrix", "factor", "missing"),
          function(exprs, phenotype, geneSet, names=NULL, weight=NULL, nBoot=100, 
                   test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE) {
            
            res <- calcGsri(exprs, phenotype, names, weight,
                            grenander, nBoot, test, testArgs, alpha)
            parms <- list(weight=weight, nBoot=nBoot, test=test, alpha=alpha,
                          grenander=grenander, testArgs=testArgs)
            object <- new("Gsri",
                          result=res$result, pval=list(res$pval), cdf=list(res$pval),
                          parms=parms)
            
            return(object)
          })

setMethod("gsri",
          signature("ExpressionSet", "factor", "missing"),
          function(exprs, phenotype, geneSet, names=NULL, weight=NULL, nBoot=100, 
                   test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE) {

            object <- gsri(exprs(exprs), phenotype, names=names, weight=weight,
                           nBoot=nBoot, test=test, testArgs=testArgs, alpha=alpha,
                           grenander=grenander)
            
            return(object)
          })

setMethod("gsri",
          signature("matrix", "factor", "GeneSet"),
          function(exprs, phenotype, geneSet, names=NULL, weight=NULL, nBoot=100, 
                   test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE) {

            if(is.null(names))
              names <- setName(geneSet)
            id <- geneIds(geneSet)
            id <- intersect(id, rownames(exprs))
            object <- gsri(exprs[id, ], phenotype, names=names, weight=weight,
                           nBoot=nBoot, test=test, testArgs=testArgs, alpha=alpha,
                           grenander=grenander)

            return(object)
          })

setMethod("gsri",
          signature("ExpressionSet", "factor", "GeneSet"),
          function(exprs, phenotype, geneSet, names=NULL, weight=NULL, nBoot=100, 
                   test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE) {

            object <- gsri(exprs(exprs), phenotype, geneSet, names=names, weight=weight,
                           nBoot=nBoot, test=test, testArgs=testArgs, alpha=alpha,
                           grenander=grenander)
            
            return(object)
          })

setMethod("gsri",
          signature("matrix", "factor", "GeneSetCollection"),
          function(exprs, phenotype, geneSet, names=NULL, weight=NULL, nBoot=100, 
                   test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE) {

            if(is.null(names))
              names <- names(geneSet)
            res <- sapply(geneSet, gsri, exprs=exprs, phenotype=phenotype, name=NULL,
                          weight=weight, nBoot=nBoot, grenander=grenander, test=test,
                          testArgs=testArgs, alpha=alpha)

            object <- new("Gsri",
                          result=as.data.frame(do.call(rbind, lapply(res, getGsri)), row.names=names),
                          pval=sapply(res, getPval),
                          cdf=sapply(res, getCdf),
                          parms=getParms(res[[1]])
                          )

            return(object)
          })

setMethod("gsri",
          signature("ExpressionSet", "factor", "GeneSetCollection"),
          function(exprs, phenotype, geneSet, names=NULL, weight=NULL, nBoot=100, 
                   test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE) {

            object <- gsri(exprs(exprs), phenotype, geneSet, names=names, weight=weight,
                         nBoot=nBoot, test=test, testArgs=testArgs, alpha=alpha,
                         grenander=grenander)
            
            return(object)
          })


## getGsri ##
setGeneric("getGsri",
           function(object)
           standardGeneric("getGsri"))

setMethod("getGsri",
          signature("Gsri"),
          function(object) {
            return(object@result)
          })

## getPval ##
setGeneric("getPval",
           function(object)
           standardGeneric("getPval"))

setMethod("getPval",
          signature("Gsri"),
          function(object) {
            return(object@pval)
          })

## getCdf ##
setGeneric("getCdf",
           function(object)
           standardGeneric("getCdf"))

setMethod("getCdf",
          signature("Gsri"),
          function(object) {
            return(object@cdf)
          })

## getParms ##
setGeneric("getParms",
           function(object)
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
setMethod("summary",
          signature("Gsri"),
          function(object, ...) {
            show(object)
          })


## plot ##
setMethod("plot",
          signature("Gsri", "integer"),
          function(x, y, digits=2, ...) {
            ##ind <- which(getGsr(x)  getGsri(x)[y, ])
            if(length(y) != 1 || y > nrow(getGsri(x)))
              stop("'y' must have exactly one match in 'x'.")
            result <- getGsri(x)[y, ]
            p <- as.numeric(result[ ,1])  ## do with names
            g <- as.numeric(result[ ,4])  ## do with names
            
            graphics::plot(c(0, 1), c(p, p), col="red", type="l", 
                           lty=2, xlab="p-values", ylab="CDF(p)", main=rownames(result), 
                           xlim=c(0, 1), ylim=c(0, 1), ...)
            graphics::lines(c(0, 1), c(g, g), col="blue", lty=2)
            graphics::lines(c(0,1), p+(1-p)*c(0,1), col="gray")
            graphics::points(x@pval[[y]], x@cdf[[y]], pch=20)
            format <- paste("%s=%.", digits, "f", sep="")
            graphics::text(1, p+0.01, sprintf(format, "%RegGene", p),
                           cex=0.8, adj=c(1, 0))
            graphics::text(1, g-0.01, sprintf(format, "GSRI", g),
                           cex=0.8, adj=c(1,1))
            
          })

setMethod("plot",
          signature("Gsri", "missing"),
          function(x, y, ...) {

            if(nrow(getGsri(x)) != 1)
              stop("'x' contains more than one gene set.")
            plot(x, 1L, ...)

          })

setMethod("plot",
          signature("Gsri", "character"),
          function(x, y, ...) {

            if(length(y) != 1)
              stop("'y' refers to more than one gene set.")
            names <- rownames(getGsri(x))
            ind <- match(y, names)
            if(is.na(ind))
              stop("'y' refers to no gene set.")
            plot(x, ind, ...)
            
          })

setMethod("plot",
          signature("Gsri", "numeric"),
          function(x, y, ...) {

            plot(x, as.integer(y), ...)

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
            phenotype <- factor(unlist(strsplit(clsCont[[3]], " ")))
            if(length(phenotype) != header[1] || nlevels(phenotype) != header[2])
              warning(sprintf("%s '%s' %s", "Data in file", basename(file),
                              "is not consistent."))
            
            return(phenotype)
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
            data <- utils::read.table(file, header=TRUE, skip=2, as.is=TRUE,
                                      row.names=1, sep="\t", quote="",
                                      na.strings=c("na", ""))
            data <- as.matrix(data[ ,-1])
            if(any(dim(data) != extend))
              warning(sprintf("%s '%s' %s", "Data in file", basename(file),
                              "is not consistent."))
            
            return(data)
          })
