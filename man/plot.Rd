\name{plot-methods}
\docType{methods}

\alias{plot}
\alias{plot,Gsri,ANY-method}
\alias{plot,Gsri,integer-method}

\title{
  Plot the GSRI results
}

\description{
  Plot the p-value distribution of the GSRI estimation.
}

\usage{
plot(x, y, ...)
}

\arguments{
  \item{x}{An object of class 'Gsri' which results are plotted.}
  \item{y}{An integer or a character specifying which gene set is
    plotted. An integer is interpreted as the index of the gene set,
    while a character is interpreted as the name of the gene set. Only
    useful if several gene sets were used in the analysis before.}
  \item{...}{Optional arguments used for plotting.}
}

\details{
  The method plots the distribution of p-values along with the final fit
  for the uniform component of the distribution.
}

\author{
  Kilian Bartholome, Julian Gehring
}

\seealso{
  \code{\link{gsri}}
}

\examples{
\dontrun{
plot(object)
}
}

\keyword{utilities}
\keyword{methods}
\keyword{hplot}
