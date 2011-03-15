\name{sortGsri}
\docType{methods}

\alias{sortGsri}
\alias{sortGsri-methods}
\alias{sortGsri,Gsri-method}

\title{
  Methods for the Gene Set Regulation Index (GSRI)
}

\description{
  Sort the result.
}

\usage{
sortGsri(x, names, decreasing=TRUE, na.last=NA, ...)
}

\arguments{
  \item{x}{Object of class 'Gsri', as obtained from the 'gsri' methods,
    whose results should be sorted.}
  \item{names}{Name(s) of the column(s) along which the results of 'x'
    should be sorted. If a vector the results are sorted first along the
    first name, then the second, and so on.}
  \item{decreasing}{Logical indicating whether the sorting should be in
    decreasing (default) or ascending order.}
  \item{na.last}{How 'NA' values in the results should be treated. See 'sort'}
  \item{...}{Additional arguments}
}

\section{Methods}{
  \describe{
    \item{
      \code{signature(x="Gsri", names="ANY")}
    }{bla}
  }
}

\details{
  This function calculates the number of differentially expressed genes
  for a single gene set, with exprs and gene set taken from the
  workspace.
  
  From bootstrapping the group samples the (1-alpha)*100\% quantile of
  the distribution of the estimated number of differentially expressed
  genes is obtained. The Gene Set Regulation Index (GSRI) is defined as
  the 5\% quantile and indicates, that with a probability of 95\% more
  than GSRI genes in the gene set are differentially expressed.

  This index can be employed to test the hypothesis whether at least one
  gene in a set is regulated and to compare and rank the regulation of
  different gene sets.
}

\value{
  An object of class 'Gsri' with the slots:

  \item{result}{Result of the estimation.}
  \item{pval}{List with the p-values for each gene set.}
  \item{cdf}{List with the cumulative denstiy for each gene set.}
  \item{parms}{Parameter values used in the analysis.}
}

\author{
  Kilian Bartholome, Julian Gehring
}

\note{
  Usage of the Grenander estimate is based on the assumption about the
  concave shape of the cumulative density distribution. It reduces the
  variance, i.e.  makes the approach more stable especially for small
  gene sets. On the downside the number of significant genes is
  overestimated for few regulated genes. A conservative solution of this
  trade-off would be to rank the gene-sets with and without Grenander
  estimate and to choose the lowest rank for each gene-set.  Please note
  that the Grenander estimate does not allow duplicates in the p-values.
  If this occurs in a data set, an error message will be printed and the
  analysis should be repeated without the Grenander estimate.

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

\seealso{
  \code{\link{GSRI}}
}

\examples{
data <- matrix(rnorm(100*20), 100)
groups <- factor(rep(1:2, each=10))

\dontrun{
## Simulate expression data for a gene set of:
## 100 genes, 40 samples (20 treatment, 20 control)
## and 30 regulated genes
expdata <- matrix(rnorm(4000,mean=0),nrow=100,ncol=40)
expdata[1:30,1:20] <- rnorm(600,mean=1)
data <- data.frame(expdata)
groups <- c(rep(0,20),rep(1,20))
geneSetName <- "Test Gene Set"

## Estimate the number of differentially expressed genes
res <- gsri(data, groups, geneSetName)
}
}
\keyword{distribution}
\keyword{htest}
