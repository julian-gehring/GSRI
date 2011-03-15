\name{read-functions}
\docType{methods}

\alias{readCls}
\alias{readGct}
\alias{readCls,character-method}
\alias{readGct,character-method}

\title{
  Import of '.cls' and '.gct' files
}

\description{
  Import the groups from '.cls' files and the expression data from
  '.gct' files.
}

\usage{
readCls(file, ...)
readGct(file, ...)
}

\arguments{
  \item{file}{Character vector specifying the path of the file to be
    read in.}
  \item{...}{Optinal arguments}
}

\details{
  With these methods the expression data and the assignment of the
  samples to groups can be read from '.cls' (categorical class) and
  '.gct' (gene cluster text) files, respectively. Details on the
  specific formats can be found at
  \url{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}.

  Please note that the 'readCls' method reads only categorical class
  labels, not continuous ones.
}

\value{
  For a '.cls' file, a factor containing the groups.

  For a '.gct' file, a matrix containing the expression intensities, with
  rows corresponding to genes and columns to samples.
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
exprs <- readGct(file)

groups <- readCls(file)
}
}

\keyword{io}