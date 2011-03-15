rowt <- function(exprs, groups, testArgs)
  genefilter::rowttests(exprs, groups)$p.value


rowF <- function(exprs, groups, testArgs)
  genefilter::rowFtests(exprs, groups)$p.value
