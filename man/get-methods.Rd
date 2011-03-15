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
  Get methods for \code{Gsri} class
}

\description{
  Access and extract results from an object of class \code{Gsri}.
}

\usage{
getGsri(object, index, ...)
getCdf(object, index, ...)
getParms(object, ...)
}

\arguments{
  \item{object}{An object of class \code{Gsri} whose results to access.}
  \item{index}{Optional argument used to subset the results of
    individual gene sets. For details, see the \code{plot} method.}
  \item{...}{Additional arguments, currently not used.}
}

\details{
  \code{getGsri} returns a data frame with the results of the GSRI analysis,
   equivalent to the display of the \code{print} and \code{show} method. For each
   gene set it contains the estimated fraction and total number of
   regulated genes, the standard deviation of the estimated fraction,
   the GSRI, and the total number of genes in the analysis.
   
   \code{getCdf} returns a data frame with ECDF of the p-values for a gene
   set.

   \code{getParms} returns a list with the parameters used in the analysis,
   including the weights (\code{weight}) for each gene in the gene set, the
   number (\code{nBoot}) of bootstraps for the calculation of the GSRI, the
   function (\code{test}) and its arguments (\code{testArgs}) used for assessing
   the differential effect between the groups, the confidence level for
   the GSRI, and the application of the Grenander estimatior
   (\code{grenander}) in the calculation of the ECDF.
}

\value{
  \code{getGsri} and \code{getCdf} return data frame with the results from the
  GSRI analysis and the ECDF, respectively. \code{getParms} returns a list
  with the parameters used in the analysis.

  If more than one gene set is accessed, a list with one element per
  gene set is returned.
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
getGsri(object)
getCdf(object)
getParms(object)
}
}

\keyword{io}
