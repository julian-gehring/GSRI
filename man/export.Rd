\name{export}
\docType{methods}

\alias{export}
\alias{export,Gsri,character-method}

\title{
  Export results to file
}

\description{
  Export the results of a \code{Gsri} object to a file.
}

\usage{
export(object, file, ...)
}

\arguments{
  \item{object}{An object of class \code{Gsri} whose results to export.}
  \item{file}{Character vector specifying the path of the file to be
    written.}
  \item{...}{Additional arguments, currently not used.}
}

\details{
  The \code{export} method writes the results of the GSRI analysis, as
  obtained with \code{getGsri}, to a text file. The file is formatted
  with the standard parameters of the \code{write.table} function, while
  \dQuote{\\t} is used as field seperator.
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

\examples{
\dontrun{
export(gsri, file)
}
}

\keyword{io}