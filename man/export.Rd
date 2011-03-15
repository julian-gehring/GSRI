\name{export}
\docType{methods}

\alias{export}
\alias{export,Gsri,character-method}

\title{
  Export of results to file
}

\description{
  Export the results of a Gsri object to a file.
}

\usage{
export(object, file, ...)
}

\arguments{
  \item{object}{An object of class 'Gsri' which results are exported.}
  \item{file}{Character vector specifying the path of the file to be
    written.}
  \item{...}{Optinal arguments}
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
export(gsri, file)
}
}

\keyword{io}