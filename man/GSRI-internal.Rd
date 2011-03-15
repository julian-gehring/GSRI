\name{GSRI-internal}

\title{Internal functions}

\alias{GSRI-internal}
\alias{gsriBoot}
\alias{multiStat}
\alias{bootInGroups}
\alias{calcGsri}
\alias{getArgs}

\description{Internal functions of the GSRI package}

\usage{
calcGsri(exprs, groups, name, id, weights, grenander=TRUE, nBoot=100,
test=NULL, testArgs=NULL, alpha=0.05, ...)
multiStat(exprs, groups, id, index, test, testArgs)
gsriBoot(exprs, groups, weights, id, grenander, test, testArgs, nSamples)
bootInGroups(nSamples)
getArgs(name, first=NULL, last=NULL, ...)
}

\author{
  Julian Gehring

  Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}

\details{
  Internal functions of the GSRI package. Users should not call them
  directly, but rather use the 'gsri' methods.
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

\keyword{internal}
