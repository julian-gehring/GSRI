\name{get-methods}
\docType{methods}

\alias{getGsri}
\alias{getCdf}
\alias{getParms}
\alias{getGsri,Gsri,ANY-method}
\alias{getGsri,Gsri,missing-method}
\alias{getCdf,Gsri,ANY-method}
\alias{getCdf,Gsri,missing-method}
\alias{getParms,Gsri-method}

\title{
  Export of results to file
}

\description{
  Export the results of a Gsri object to a file.
}

\usage{
getGsri(object, index, ...)
getCdf(object, index, ...)
getParms(object, ...)
}

\arguments{
  \item{object}{An object of class 'Gsri' which results are exported.}
  \item{index}{Optional argument used to subset the results.}
  \item{...}{Additional arguments, currently not used.}
}

\details{
  The method writes a table with results to a text file. The output is
  formatted as given by the standard parameters of 'write.table' and can
  be read in directly by 'read.table' without further specification
  about the format.
}

\author{
  Kilian Bartholome, Julian Gehring
}

\seealso{
  \code{\link{gsri}}
}

\examples{
\dontrun{
getGsri(object)
getCdf(object)
getParms(object)
}
}

\keyword{io}
