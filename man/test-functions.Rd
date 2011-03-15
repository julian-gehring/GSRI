\name{test-functions}

\alias{rowt}
\alias{rowF}
\alias{limmat}

\title{
  Statistical test functions
}

\description{
  Estimate the effect in gene expression between groups of microarray
  replicates.
}

\usage{
rowt(exprs, groups, id, index, testArgs)
rowF(exprs, groups, id, index, testArgs)
limmat(exprs, groups, id, index, testArgs)
}

\arguments{
  \item{exprs}{A matrix of expression values of size n x m, with rows
    representing the genes and columns representing the samples. The
    structure is the same as for the 'exprs' argument of the 'gsri'
    method.}
  \item{groups}{A factor with the length m, specifying the groups
    of the corresponding samples in 'exprs'.  The structure is the same
    as for the 'exprs' argument of the 'groups'  method.}
  \item{id}{Index vector for the rows of 'exprs' which are part of the
    current gene set.}
  \item{index}{Index for the columns of 'exprs', such that 'exprs[
    ,index] yields the bootstrapped expression matrix. Similar to the
    'index' arguments from 'boot' of the 'boot' package.}
  \item{testArgs}{Optional arguments passed to the test function.}
}

\details{
  With the t-test and the F-test, two widely used statistical tests are
  available in this package. To allow a fast computation the
  implementations from the \pkg{genefilter} package is used.

  It is also possible to use custom test statistics for assessing the
  differential effect between groups for each gene. In this case the
  function is passed as the 'test' argument to the 'gsri' method, while
  additional parameters for the function can be passed via the 'testArgs'
  argument. The defined function is required to be called as
  
  \code{function(exprs, groups, id, index, testArgs)},

  with 'exprs' the matrix of expression intensities of the microarray
  and 'groups' the factor of group labels, with the same structure as
  those passed initially to the 'gsri' method. The vector 'id' contains
  the indices of the genes part of the current gene set and is used to
  subset the expression intensities if necessary. The function has to
  return one p-value for each gene in the gene set indicating its
  differential effect. The vector 'index' contains the indicies of the
  samples for the bootstrapping. Applying 'index' on the expression
  matrix in the form of \code{exprs[ ,index]} generates the bootstrapped
  data set.
  
  For details on how to define and use your custom test functions,
  please refer to the 'examples' section or the vignette of this
  package.
}

\value{
  A vector of p-values, indicating the significance of the differential
  effect between groups for each gene.
}

\author{
  Julian Gehring

  Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}

\seealso{
  Package:
  \code{\link[GSRI]{GSRI-package}}
  
  Class:
  \code{\linkS4class{Gsri}}
  
  Methods:
  \code{\link[GSRI]{gsri}}  
  \code{\link[GSRI]{getGsri}}
  \code{\link[GSRI]{getCdf}}
  \code{\link[GSRI]{getParms}}
  \code{\link[GSRI]{export}}
  \code{\link[GSRI]{sortGsri}}
  \code{\link[GSRI]{plot}}
  \code{\link[GSRI]{show}}
  \code{\link[GSRI]{summary}}
  \code{\link[GSRI]{readCls}}
  \code{\link[GSRI]{readGct}}
}

\examples{
## A simple example for a custom test function using a linear model.
## Note that for two groups this is equivalent to a t-test with equal variances.


testFcn <- function(exprs, groups, id, index, testArgs) {

  stat <- function(e, g) {
  m <- lm(e ~ g)
  pval <- summary(m)$coefficients[2,4]
  }

pvals <- apply(exprs[id,index], 1, statFcn, groups)
return(pvals)
}

\dontrun{
g <- gsri(exprs, groups, test=testFcn)
}
}

\keyword{htest}