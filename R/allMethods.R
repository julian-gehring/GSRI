setGeneric("gsri",
           function(exprs, phenotype, geneSet, ...)
           standardGeneric("gsri"))

setMethod("gsri",
          signature("matrix", "factor", "missing"),
          function(exprs, phenotype, geneSet, names=NULL, weight=NULL, nBoot=100, 
                   test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE) {
            
            res <- calcGsri(exprs, phenotype, names, weight,
                            grenander, nBoot, test, testArgs, alpha)
            ##object <- new("Gsri", result=res$result, pval=list(res$pval))
            parms <- list(weight=weight, nBoot=nBoot, test=test, alpha=alpha,
                          grenander=grenander, testArgs=testArgs)
            object <- new("Gsri", result=res, parms=parms)
            
            return(object)
          })

setMethod("gsri",
          signature("ExpressionSet", "factor", "missing"),
          function(exprs, phenotype, geneSet, names=NULL, weight=NULL, nBoot=100, 
                   test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE) {

            res <- calcGsri(exprs(exprs), phenotype, names, weight,
                            grenander, nBoot, test, testArgs, alpha)
            object <- new("Gsri", result=res)
            
            return(object)
          })

setMethod("gsri",
          signature("matrix", "factor", "GeneSet"),
          function(exprs, phenotype, geneSet, names=NULL, weight=NULL, nBoot=100, 
                   test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE) {
            
            id <- geneIds(geneSet)
            id <- intersect(id, rownames(exprs))
            res <- calcGsri(exprs[id, ], phenotype, names, weight,
                            grenander, nBoot, test, testArgs, alpha)            
            object <- new("Gsri", result=res)
            return(object)
          })

setMethod("gsri",
          signature("ExpressionSet", "factor", "GeneSet"),
          function(exprs, phenotype, geneSet, names=NULL, weight=NULL, nBoot=100, 
                   test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE) {
            
            id <- geneIds(geneSet)
            res <- calcGsri(exprs(exprs[id, ]), phenotype, names, weight,
                            grenander, nBoot, test, testArgs, alpha)            
            object <- new("Gsri", result=res)
            return(object)
          })

setMethod("gsri",
          signature("ANY", "factor", "GeneSetCollection"),
          function(exprs, phenotype, geneSet, names=NULL, weight=NULL, nBoot=100, 
                   test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE) {
            
            extract <- function(geneSet, data, phenotype, name, weight, grenander,
                                nBoot, test, testArgs, alpha) {
              id <- geneIds(geneSet)
              res <- calcGsri(exprs[id, ], phenotype, names, weight,
                              grenander, nBoot, test, testArgs, alpha)
            }
            res <- sapply(geneSet, extract, data=exprs, phenotype=phenotype, name=names,
                          weight=weight, nBoot=nBoot, grenander=grenander, test=test,
                          testArgs=testArgs, alpha=alpha)
            res <- as.data.frame(t(res))
            if(is.null(names))
              names <- names(gsc)
            rownames(res) <- names
            object <- new("Gsri", result=res)
            
            return(object)
          })


setMethod("show",
          signature("Gsri"),
          function(object) {
            print(object@result)
          })

setGeneric("getGsri",
           function(object)
           standardGeneric("getGsri"))

setMethod("getGsri",
          signature("Gsri"),
          function(object) {
            object@result
          })


setGeneric("getParms",
           function(object)
           standardGeneric("getParms"))

setMethod("getParms",
          signature("Gsri"),
          function(object) {
            data.frame(object@parms, row.names="")
          })

setMethod("plot",
          signature("Gsri", "ANY"),
          function(x, y, digits=2, ...) {
            ##ind <- which(getGsr(x)  getGsri(x)[y, ])
            ind <- y
            if(length(ind) != 1)
              stop("'y' must have exactly one match in 'x'.")
            result <- getGsri(x)[ind, ]
            p <- as.numeric(result[1])
            g <- as.numeric(result[3])
            
            graphics::plot(c(0, 1), c(p, p), col="red", type="l", 
                           lty=2, xlab="p-values", ylab="CDF(p)", main=rownames(result), 
                           xlim=c(0, 1), ylim=c(0, 1), ...)
            graphics::lines(c(0, 1), c(g, g), col="blue", lty=2)
            graphics::lines(c(0,1), p+(1-p)*c(0,1), col="gray")
            graphics::points(x@pval[[ind]], x@cdf[[ind]], pch=20)
            format <- paste("%s=%.", digits, "f", sep="")
            graphics::text(1, p+0.01, sprintf(format, "%RegGene", p),
                           cex=0.8, adj=c(1, 0))
            graphics::text(1, g-0.01, sprintf(format, "GSRI", g),
                           cex=0.8, adj=c(1,1))
            
          })


setGeneric("export",
           function(object, file, ...)
           standardGeneric("export"))

setMethod("export",
          signature("Gsri", "character"),
          function(object, file, digits=Inf) {
            result <- round(getGsri(object), digits)
            write.table(result, file, sep="\t")
          })


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
