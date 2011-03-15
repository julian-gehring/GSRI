\name{read-functions}
\docType{methods}

\alias{readCls}
\alias{readGct}
\alias{readCls,character-method}
\alias{readGct,character-method}

\title{
  Import of 'cls' and 'gct' files
}

\description{
  Import the groups from 'cls' files and the expression data from
  'gct' files.
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
  The input files should have the format typical for those file
  types. Details on the specific formats can be found at
  \url{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}.
}

\value{
  For a cls file, a factor containing the groups.

  For a gct file, a matrix containing the expression intensities.
}

\author{
  Kilian Bartholome, Julian Gehring
}

\seealso{
  \code{\link{gsri}}
}

\examples{
\dontrun{
exprs <- readGct(file)

groups <- readCls(file)
}
}

\keyword{io}