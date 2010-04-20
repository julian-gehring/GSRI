\name{gsri-internal}
\title{Internal GSRI functions}
\alias{gsri-internal}
\alias{getCDF}
\alias{getPvalues}
\alias{grenanderInterp}
\alias{gsriBoot}
\alias{loadCls}
\alias{loadGct}
\alias{multiStat}
\alias{slopeFast}

\description{Internal GSRI functions}

\usage{
getCDF(pvals, useGrenander)
getPvalues(data, d, phenotype, test, testArgs)
grenanderInterp(x, y)
gsriBoot(inputData, d, phenotype, useGrenander, test, testArgs)
loadCls(fileName)
loadGct(fileName)
multiStat(data, phenotype, test, testArgs)
slopeFast(x, y)
}
\author{Kilian Bartholomé, Julian Gehring}

\details{
  Internal functions of the GSRI package.
  These are not to be called by the user.
}

\references{
Functions of the following packages are used in this package:

Korbinian Strimmer (2009).
fdrtool: Estimation and Control of (Local) False Discovery Rates.
R package version 1.2.6.
\url{http://CRAN.R-project.org/package=fdrtool}

R. Gentleman, V. Carey, W. Huber and F. Hahne.
genefilter: methods for filtering genes from microarray experiments.
R package version 1.28.2.
}


\seealso{
\code{\link{gsri}}

\code{\link{gsriFromFile}}

\code{\link{GSRI}}
}
\keyword{internal}

