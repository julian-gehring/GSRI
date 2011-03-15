\name{export}
\docType{methods}

\alias{export}
\alias{export,Gsri,character-method}

\title{
  Export results to file
}

\description{
  Export the results of a 'Gsri' object to a file.
}

\usage{
export(object, file, ...)
}

\arguments{
  \item{object}{An object of class 'Gsri' whose results to export.}
  \item{file}{Character vector specifying the path of the file to be
    written.}
  \item{...}{Optional arguments}
}

\details{
  The 'export' method writes the results of the GSRI analysis, as
  obtained with 'getGsri', to a text file. The file is formatted with
  the standard parameters of the 'write.table' function, while '\\t' is
  used as field seperator.
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