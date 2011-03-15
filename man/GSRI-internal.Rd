\name{gsri-internal}

\title{Internal GSRI functions}

\alias{gsri-internal}
\alias{getPvalues}
\alias{gsriBoot}
\alias{multiStat}
\alias{bootWithinGroups}
\alias{calcGsri}

\description{Internal GSRI functions}

\usage{
getPvalues(data, d, phenotype, test, testArgs)
gsriBoot(data, index, phenotype, cweight, grenander, test, testArgs, ...)
multiStat(data, label, weight, grenander, test, ...)
bootWithinGroup(exprs, label)
calcGsri(data, phenotype, name, weight, grenander=TRUE, nBoot=100, test="ttest", testArgs=NULL, alpha=0.05, ...)
}

\author{
  Kilian Bartholome, Julian Gehring
}

\details{
  Internal functions of the GSRI package.
  These are not to be called by the user.
}

\references{
  Functions of the following packages are used in this package:

  Korbinian Strimmer (2009).
  fdrtool: Estimation and Control of (Local) False Discovery Rates.
  \url{http://CRAN.R-project.org/package=fdrtool}

  R. Gentleman, V. Carey, W. Huber and F. Hahne.
  genefilter: methods for filtering genes from microarray experiments.
}

\seealso{
  \code{\link{gsri}}

  \code{\link{GSRI}}
}

\keyword{internal}
