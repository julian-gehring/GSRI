rowt <- function(exprs, groups, id, testArgs)
  genefilter::rowttests(exprs[id, ], groups)$p.value


rowF <- function(exprs, groups, id, testArgs)
  genefilter::rowFtests(exprs[id, ], groups)$p.value
