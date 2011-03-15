setClass("Gsri",
         representation(result="data.frame",
                        pval="list",
                        cdf="list",
                        parms="list"))

setValidity("Gsri",
            function(object) {
              res <- TRUE
              if(length(object@pval) != length(object@cdf)) {
                res <- "pval and cdf do not have same length"
                return(res)
              }
              if(nrow(object@result) != length(object@pval)) {
                res <- "pval and result do not have same length"
                return(res)
              }
              pval <- unlist(object@pval)
              cdf <- unlist(object@cdf)
              if(length(object@pval) > 0 || length(object@cdf) > 0) {
                if(min(pval, na.rm=TRUE) < 0 || max(pval, na.rm=TRUE) > 1) {
                  res <- "pval not valid"
                  return(res)
                }
                if(min(cdf, na.rm=TRUE) < 0  || max(cdf, na.rm=TRUE) > 1) {
                  res <- "cdf not valid"
                  return(res)
                }
              }
              return(res)
            })
