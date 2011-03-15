\name{GSRI-internal}

\title{Internal functions}

\alias{GSRI-internal}
\alias{gsriBoot}
\alias{multiStat}
\alias{bootInGroups}
\alias{calcGsri}
\alias{getArgs}

\description{Internal functions of the \pkg{GSRI} package}

\usage{
calcGsri(exprs, groups, name, id, weights, grenander=TRUE, nBoot=100,
test=NULL, testArgs=NULL, alpha=0.05, verbose=FALSE, ...)
multiStat(exprs, groups, id, index, test, testArgs)
gsriBoot(exprs, groups, weights, id, grenander, test, testArgs, nSamples)
bootInGroups(nSamples)
getArgs(name, first=NULL, last=NULL, ...)
}

\details{
  Internal functions of the \pkg{GSRI} package. Users should not call them
  directly, but rather use the \code{gsri} methods.
}

\references{
  The \pkg{GSRI} package uses the functionality of the following
  packages:

  Julian Gehring, Clemens Kreutz, Jens Timmer:
  les: Identifying Loci of Enhanced Significance in Tiling Microarray Data
  \url{http://bioconductor.org/help/bioc-views/release/bioc/html/les.html}

  Korbinian Strimmer:
  fdrtool: Estimation and Control of (Local) False Discovery Rates.
  \url{http://CRAN.R-project.org/package=fdrtool}

  Robert Gentleman, Vincent J. Carey, Wolfgang Huber, Florian Hahne:
  genefilter: methods for filtering genes from microarray experiments.
  \url{http://bioconductor.org/help/bioc-views/release/bioc/html/genefilter.html}
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

\keyword{internal}
