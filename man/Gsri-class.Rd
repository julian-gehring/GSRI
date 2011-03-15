\name{Gsri-class}

\docType{class}

\alias{Gsri}
\alias{Gsri-class}
\alias{show}
\alias{show,Gsri-method}
\alias{summary}
\alias{summary,Gsri-method}
\title{
  Class \code{Gsri}
}

\description{
  Objects of the class \code{Gsri} contain the results of the GSRI analysis.
}

\section{Objects from the class}{
  Objects of class \code{Gsri} are returned by the \code{gsri} methods.
}

\section{Slots}{
  \describe{
    \item{\code{result}:}{Data frame containing the results of the GSRI
      estimation, with one row for each gene set and the columns:
      \describe{
	\item{\code{pRegGenes}:}{Fraction of regulated genes in the gene set}
	\item{\code{pRegGenesSd}:}{Standard deviation of \code{pRegGenes} obtained
	  from bootstrapping.}
	\item{\code{nRegGenes}:}{Total number of regulated genes in the gene
	  set.}
	\item{\code{GSRI('alpha'\%)}:}{Gene Set Regulation Index, corresponding to
	  the \sQuote{alpha}\% quantile of the bootstrapped distribution.}
	\item{\code{nGenes}:}{Total number of genes in the gene set.}
      }
    }
    \item{\code{cdf}:}{List of data frames containing the ECDF of the
      p-values. Each data frame covers one gene set, with the columns:
      \describe{
	\item{\code{pval}:}{P-values obtained from the \code{test} function.}
	\item{\code{cdf}:}{Empirical cumulative density.}
      }
    }
    \item{\code{parms}:}{List containing the parameter values used in
      the analysis, with the elements:
      \describe{
	\item{\code{weight}:}{Weights for each gene in the gene set}
	\item{\code{nBoot}:}{Number of bootstraps for the calculation of the GSRI}
	\item{\code{test}:}{Statistical test function}
	\item{\code{alpha}:}{Confidence level for the GSRI}
	\item{\code{grenander}:}{Application of the Grenander estimatior in the
	  calculation of the ECDF}
	\item{\code{testArgs}:}{Optional arguments for \code{test} function}
      }
    }
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
      
      Assess the degree of differential effect in the expression data.}
  }

  
  Visualization:
  \describe{
    \item{plot:}{
      \code{signature(x="Gsri", y=ANY)}
    }
  }
  Plot the empirical density of p-values and the corresponding estimated effect.

  
  Export to file:
  \describe{
    \item{export:}{
      \code{signature(object="Gsri", file="character")}
    }
  }

  
  Get methods:
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

  
  Show:
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
showClass("Gsri")
}

\keyword{classes}
\keyword{methods}
