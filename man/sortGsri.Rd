\name{sortGsri}
\docType{methods}

\alias{sortGsri}
\alias{sortGsri-methods}
\alias{sortGsri,Gsri-method}

\title{
  Sort GSRI results
}

\description{
  Sort the results of an 'Gsri' object.
}

\usage{
sortGsri(x, names, decreasing=TRUE, na.last=NA, ...)
}

\arguments{
  \item{x}{Object of class 'Gsri' whose results to sort.}
  \item{names}{Columns along which the results of 'x' should be sorted,
    eighter a character vector with the names of the columns or an
    integer vector with the index of the columns. If the vector has
    several elements, sorting is performed along all of them, starting
    with the first and using subsequent ones to break existing ties. If
    names is not specified the results are sorted according to 'pRegGenes'.}
  \item{decreasing}{Logical indicating whether the sorting should be in
    decreasing (default) or ascending order, see 'sort'.}
  \item{na.last}{How 'NA' values in the results should be treated, see 'sort'.}
  \item{...}{Additional arguments, currently not used.}
}

\section{Methods}{
  \describe{
    \item{
      \code{signature(x="Gsri", names="ANY")}
    }{}
  }
}

\details{
  nothing
}

\value{
  An object of class 'Gsri', with sorted slots 'result' and 'cdf'.
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
sortGsri(object, c("pRegGenes", "nGenes"))
sortGsri(object, c(1, 5))
}
}

\keyword{utils}
