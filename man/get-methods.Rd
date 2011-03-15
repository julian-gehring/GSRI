\name{get-methods}
\docType{methods}

\alias{getGsri}
\alias{getParms}
\alias{getGsri,Gsri-method}
\alias{getParms,Gsri-method}

\title{
  Export of results to file
}

\description{
  Export the results of a Gsri object to a file.
}

\usage{
getGsri(object)
getParms(object)
}

\arguments{
  \item{object}{An object of class 'Gsri' which results are exported.}
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
getParms(object)
}
}

\keyword{io}
