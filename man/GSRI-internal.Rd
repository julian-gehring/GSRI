\name{GSRI-internal}

\title{Internal GSRI functions}

\alias{GSRI-internal}
\alias{gsriBoot}
\alias{multiStat}
\alias{bootWithinGroups}
\alias{calcGsri}

\description{Internal GSRI functions}

\usage{
gsriBoot(exprs, groups, weights, id, grenander, test, testArgs)
multiStat(exprs, groups, id, grenander, test, testArgs)
bootWithinGroup(exprs, groups, weights)
calcGsri(exprs, groups, name, id, weights, grenander=TRUE, nBoot=100, test="ttest", testArgs=NULL, alpha=0.05, ...)
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
