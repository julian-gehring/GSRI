rowt <- function(data, phenotype, ...)
  genefilter::rowttests(data, phenotype)$p.value


rowF <- function(data, phenotype, ...)
  genefilter::rowFtests(data, phenotype)$p.value
