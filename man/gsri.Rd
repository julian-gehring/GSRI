\name{gsri}

\docType{methods}

\alias{gsri}
\alias{gsri-methods}
\alias{gsri,ExpressionSet,factor,GeneSet-method}
\alias{gsri,ExpressionSet,factor,GeneSetCollection-method}
\alias{gsri,ExpressionSet,factor,missing-method}
\alias{gsri,matrix,factor,GeneSet-method}
\alias{gsri,matrix,factor,GeneSetCollection-method}
\alias{gsri,matrix,factor,missing-method}

\title{
  Methods for the Gene Set Regulation Index (GSRI)
}

\description{
  Estimate the number of differentially expressed genes in gene sets.
}

\usage{
gsri(exprs, groups, geneSet, names=NULL, weight=NULL, nBoot=100,
test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE, ...)
}

\arguments{
  \item{exprs}{Matrix or object of class \code{ExpressionSet} containing
    the expression intensities of the microarray. If a matrix the rows
    represent the genes and the columns the samples, respectively.}
  \item{groups}{Factor with the assignments of the microarray samples to
    the groups, along which the differential effect should be
    estimated. Must have as many elements as \code{exprs} has samples.}
  \item{geneSet}{Optional object of class \code{GeneSet} or
    \code{GeneSetCollection} defining the gene set(s) used for the
    analysis. If missing all genes of \code{exprs} are considered to be part
    of the gene set. If an object of class \code{GeneSet} only these genes
    are considered to be part of the gene set. If an object of class
    \code{GeneSetCollecton} the analysis is performed for each gene set of
    the collection individually.}
  \item{names}{Optional character vector with the names of the gene
    set(s). If missing the names are taken from the \code{geneSet}
    argument. Has to have as many unique elements as gene sets in the
    analysis.}
  \item{weight}{Optional numerical vector of weights specifying the
    certainty a gene is part of the gene set. If \code{NULL} all genes are
    assumed to have the same weight. Please note that the weights are
    defined in a relative way and thus any kind of positive
    weights is feasable. Must have as many elements as eighter
    the genes defined in \code{geneSet} or in \code{exprs}.}
  \item{nBoot}{Integer with the number of bootstrap samples to be drawn
    in the calculation of the GSRI (default: 100).}
  \item{test}{A function defining the statistical test used to assess
    the differential effect between the groups which are given by the
    \code{groups} argument. In this package, a t-test (rowt) and an F-test
    (rowF) are already supplied, with \code{rowt} being the
    default. Additionally, a custom test function can be used in order
    to be able to include any feasible statistical test in the
    analysis. For details, please see the \sQuote{details} section.}
  \item{testArgs}{List of optional arguments used by the \code{test}
    function. For details, please see the \sQuote{details} section and the help
    for \code{test-functions}.}
  \item{alpha}{Single numeric specifying the confidence level for the
    GSRI. The estimated GSRI is the lower bound of the
    (1-\eqn{\alpha}{\sQuote{alpha}})*100\% confidence interval obtained from
    bootstrapping.}
  \item{grenander}{Logical about whether the modified Grenander
    estimator for the cumulative density should be used instead of a
    centered ECDF. By default the modified Grenander estimator is
    used. For more information, please see the \sQuote{details} section.}
  \item{...}{Additional arguments, including:
    \describe{
      \item{nCores:}{Integer setting the number of cores used for the
	computation in combination with the \pkg{multicore} package for
	a \code{GeneSetCollection}. For details, please see the \sQuote{details}
	section.}
      \item{verbose:}{Logical indicating whether the progress of the
	computation for a \code{GeneSetCollection} should be printed to
	screen (default: FALSE).}
    }
  }
}

\section{Methods}{
  Analysis for all genes of \code{exprs} part of the gene set:
  \describe{
    \item{
      \code{signature(exprs="matrix", groups="factor", geneSet="missing")}
    }{}

    \item{
      \code{signature(exprs="ExpressionSet", groups="factor", geneSet="missing")}
    }{}
  }

  Analysis for one gene set, defined as an object of class \code{GeneSet}:
  \describe{
    \item{
      \code{signature(exprs="matrix", groups="factor", geneSet="GeneSet")}
    }{}
    
    \item{
      \code{signature(exprs="ExpressionSet", groups="factor", geneSet="GeneSet")}
    }{}
  }

  Analysis for several gene sets, defined as an object of class
  \code{GeneSetCollection}:
  \describe{
    \item{
      \code{signature(exprs="matrix", groups="factor", geneSet="GeneSetCollection")}
    }{}
    
    \item{
      \code{signature(exprs="ExpressionSet", groups="factor", geneSet="GeneSetCollection")}
    }{}
  }
  In this case parallel computing capabilities provided by the \pkg{multicore}
  package may be available, depending on the platform.
}

\details{
  The \code{gsri} method estimates the degree of differential expression in
  gene sets. By assessing the part of the distribution of p-values
  consistent with the null hypothesis the number of differentially
  expressed genes is calculated.

  Through non-parametric fitting of the uniform component of the p-value
  distribution, the fraction of regulated genes \eqn{r}{\sQuote{r}} in a gene
  set is estimated. The GSRI \eqn{\eta}{\sQuote{eta}} is then defined as the
  \eqn{\alpha*100}{\sQuote{alpha*100}}\%-quantile of the distribution of
  \eqn{r}{\sQuote{r}}, obtained from bootstrapping the samples within the
  groups. The index indicates that with a probability of
  \eqn{(1-\alpha)}{(1-\sQuote{alpha})}\% more than a fraction of
  \eqn{\eta}{\sQuote{eta}} genes in the gene set is differentially
  expressed. It can also be employed to test the hypothesis whether at
  least one gene in a gene set is regulated. Further, different gene
  sets can be compared or ranked according to the estimated amount of
  regulation.

  Assessing the differential effect is based on p-values obtained from
  statistical testing at the level of individual genes between the
  groups. The GSRI approach is independent of the underlying test and
  can be chosen according to the experimental design. With the t-test
  (\code{rowt}) and F-test (\code{rowF}) two widely used statistical test are
  already part of the package. Additional tests can easily used which
  are passed with the \code{test} argument to the \code{gsri} method. For details
  on how to implement custom test functions, please refer to the help of
  \code{rowt} and \code{rowF} or the vignette of this package.

  The GSRI approach further allows weighting the influence of individual
  genes in the estimation. This can be beneficial including for example
  the certainty that genes are part of a certain gene set derived from
  experimental findings or annotations.

  Defining gene sets is available through the \pkg{GSEABase} package which
  provides the \code{GeneSet} and \code{GeneSetCollection} classes a single or
  multiple gene sets, respectively. This ensures a powerful approach for
  obtaining gene sets from data objects, data bases, and other
  bioconductor packages. For details on how to define or retrieve gene
  sets, please refer to the documentation of the \pkg{GSEABase} package,
  with a special focus on the \code{GeneSet} and \code{GeneSetCollection} classes.

  The distribution of the p-values of a gene set is assessed in the
  cumulative density. In addition to a symmetrical empirical cumulative
  density function (ECDF), the modified Grenander estimator based on the
  assumption about the concave shape of the cumulative density is
  implemented and used by default. While the modified Grenander
  estimator reduces the variance and makes the approach more stable
  especially for small gene set, it underestimates the number of
  regulated genes and thus leads to conservative estimates.

  In the case that the computation is performed for several gene sets in
  the form of a \code{GeneSetCollection} object, it can be parallelized with the
  \pkg{multicore} package. Please note that this package is not available
  on all platforms. Using its capabilities requires attaching
  \pkg{multicore} prior to the calculation and specification of the \code{nCores}
  argument. For further details, please refer to the documentation of
  the \pkg{multicore} package. This may be especially relevant in the case
  that specific seed values for the bootstrapping are of interest.
}

\value{
  An object of class \code{Gsri} with the slots:
  \describe{
    \item{\code{result}:}{Data frame containing the results of the GSRI
      estimation, with one row for each gene set.}
    \item{\code{cdf}:}{List of data frames containing the ECDF of the
      p-values. Each data frame covers one gene set.}
    \item{\code{parms}:}{List containing the parameter values used in
      the analysis, with the elements.}
  }
  For details, please see the help for the \code{\linkS4class{Gsri}}
  class.
}

\note{
  The standard deviation of the estimated number of regulated genes as
  well as the GSRI are obtained through bootstrapping. Thus, the results
  for these two parameters may differ slightly for several realizations,
  especially for small numbers of bootstraps (\code{nBoot}). Setting the seed
  of the random number generator avoids this problem and yields exactly
  the same results for several realizations.
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
## Simulate expression data for a gene set of
## 100 genes, 20 samples (10 treatment, 10 control)
## and 30 regulated genes
set.seed(1)
exprs <- matrix(rnorm(100*20), 100)
exprs[1:30,1:10] <- rnorm(30*10, mean=2)
groups <- factor(rep(1:2, each=10))

## Estimate the number of differentially expressed genes
res <- gsri(exprs, groups)
}

\keyword{distribution}
\keyword{htest}
