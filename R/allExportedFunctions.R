rowt <- function(exprs, groups, ...)
  genefilter::rowttests(exprs, groups)$p.value


rowF <- function(exprs, groups, ...)
  genefilter::rowFtests(exprs, groups)$p.value
