setClass("Gsri",
         representation(result="data.frame",
                        cdf="list",
                        parms="list"))

setValidity("Gsri",
            function(object) {
              res <- TRUE
              if(nrow(object@result) != length(object@cdf)) {
                res <- "cdf and result do not have same length"
                return(res)
              }
              return(res)
            })
