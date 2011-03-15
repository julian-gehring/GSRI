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
