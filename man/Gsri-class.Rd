\name{Gsri-class}

\docType{class}

\alias{Gsri}
\alias{Gsri-class}
\alias{show,Gsri-method}
\alias{summary,Gsri-method}
\title{
  Class 'Gsri'
}

\description{
  Objects of the class 'Gsri' containing the results of the estimation
  of the Gene Set Regulation Index.
}

\section{Objects from the class}{
  Objects of class 'Gsri' are the output of the method 'gsri'.
}

\section{Slots}{
  \describe{
    \item{\code{result}:}{Data frame containing the results of the GSRI
      estimation, including the ratio and number of regulated genes, the
      total number of genes in the gene set, and the Gene Set Regulation
      Index.}
    \item{\code{cdf}:}{Data frame containing the empirical cumulative
      density for each gene set.}
    \item{\code{parms}:}{List containing the parameter values used in
      the analysis.}
  }
}

\section{Methods}{
  Analysis:
  \describe{
    \item{gsri:}{
      \code{signature(exprs="matrix", groups="factor",
	geneSet="missing")}
      
      \code{signature(exprs="ExpressionSet", groups="factor",
	geneSet="missing")}
      
      \code{signature(exprs="matrix", groups="factor",
	geneSet="GeneSet")}
      
      \code{signature(exprs="ExpressionSet", groups="factor",
	geneSet="GeneSet")}
      
      \code{signature(exprs="matrix", groups="factor",
	geneSet="GeneSetCollection")}
      
      \code{signature(exprs="ExpressionSet", groups="factor",
	geneSet="GeneSetCollection")}
      
      do the analysis}
  }

  
  Visualization:
  \describe{
    \item{plot:}{
      \code{signature(x="Gsri", y=ANY)}
    }
  }

  
  Export:
  \describe{
    \item{export:}{
      \code{signature(object="Gsri", file="character")}
    }
  }

  
  Slot access:
  \describe{
    \item{getGsri:}{
      \code{signature(object="Gsri")}
    }

    \item{getCdf:}{
      \code{signature(object="Gsri")}
    }

    \item{getParms:}{
      \code{signature(object="Gsri")}
    }
  }

  
  Report:
  \describe{
    \item{show:}{
      \code{signature(obejct="Gsri")}
    }
    \item{summary:}{
      \code{signature(obejct="Gsri")}
    }
  }
}


\author{
  Kilian Bartholome, Julian Gehring, Clemens Kreutz, Jens Timmer
}

\seealso{
  Package:
  \code{\link[GSRI]{GSRI-package}}

  Class:
  \code{\linkS4class{Gsri}}

  Methods:
  \code{\link[GSRI]{gsri}}
  \code{\link[GSRI]{plot}}
  \code{\link[GSRI]{export}}
  \code{\link[GSRI]{getGsri}}
  \code{\link[GSRI]{getParms}}
}

\examples{
\dontrun{
showClass("Gsri")
}
}

\keyword{classes}
\keyword{methods}
