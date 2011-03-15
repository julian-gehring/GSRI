\name{test-functions}

\alias{rowt}
\alias{rowF}

\title{
  Statistical test functions
}

\description{
  Estimate the effect in gene expression between groups of microarray
  replicates.
}

\usage{
rowt(exprs, groups, id, testArgs)
rowF(exprs, groups, id, testArgs)
}

\arguments{
  \item{exprs}{A matrix of expression values of size n x m, with rows
    representing the genes and columns representing the arrays.}
  \item{groups}{A factor with the length m, specifying the groups
    of the corresponding arrays in 'exprs'.}
  \item{id}{Index for the rows of 'exprs' which are part of the
    current gene set.}
  \item{testArgs}{Optinal arguments, not used in the default functions.}
}

\details{
  With the t-test and the F-test, two widely used statistical tests are
  available in this package. To allow fast computation this package uses
  the implementations from the \pkg{genefilter} package.

  It is also possible to apply user-defined tests with this method. In
  this case the function has to be called by function(exprs, groups,
  testArgs). 'exprs' and 'groups' are of class 'matrix' and 'factor',
  respectively. 'testArgs' will only be passed to the function if it is
  defined. In general all methods that compute p-values are suitable. The
  function has to return a vector with one p-value per gene. For details
  on how to use your own test functions please refer to the vignette of
  this package.
}

\value{
  A vector of p-values, indicating the significance of the effect for
  each gene.
}

\author{
  Kilian Bartholome, Julian Gehring
}

\seealso{
  \code{\link{gsri}}
}

\examples{
statFcn <- function(d, p) {
    m <- lm(d ~ p)
    pval <- summary(m)$coefficients[2,4]
}

testFcn <- function(exprs, groups, ...) {
    pvals <- apply(exprs, 1, statFcn, groups)
    return(pvals)
}
}

\keyword{htest}