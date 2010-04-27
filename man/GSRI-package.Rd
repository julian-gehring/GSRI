\name{GSRI-package}
\alias{GSRI-package}
\alias{GSRI}
\docType{package}
\title{Gene Set Regulation Index (GSRI)}

\description{
  This package estimates the number of differentially
  expressed genes in gene sets with the Gene Set Regulation Index
  (GSRI).
}

\details{
  
  \tabular{ll}{
    Package: \tab GSRI\cr
    Type: \tab Package\cr
    Version: \tab 1.1.1\cr
    Date: \tab 2010-04-27\cr
    License: \tab GPL-2\cr
    LazyLoad: \tab yes\cr
  }
  This method calculates the number of differentially expressed genes in
  a gene set. It does not require a cut-off value for the distinction
  between regulated and unregulated genes. The approach is based on the
  fact that p-values obtained from a statistical test are uniformly
  distributed under the null hypothesis and accumulated around zero for
  the alternative hypothesis.

  The Gene Set Regulation Index (GSRI) is the 5\%-quantile of the
  distribution of the estimated number of differentially expressed genes
  obtained from bootstrapping the group samples. It indicates, that with
  a probability of 95\% more than GSRI genes in the gene set are
  differentially expressed. This index can also be employed to test the
  hypothesis whether at least one gene in a set is regulated and to
  compare and rank the regulation of different gene sets, see Bartholomé
  et al. 2009.
}

\author{
  Kilian Bartholomé, Julian Gehring

  Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}

\references{
  Kilian Bartholomé, Clemens Kreutz, and Jens Timmer:
  Estimation of gene induction enables a relevance-based ranking of gene sets,
  Journal of Computational Biology: A Journal of Computational Molecular Cell 
  Biology 16, no. 7 (July 2009): 959-967.
  \url{http://www.liebertonline.com/doi/abs/10.1089/cmb.2008.0226}

  Functions of the following packages are used in this package:

  Korbinian Strimmer (2009).
  fdrtool: Estimation and Control of (Local) False Discovery Rates.
  R package version 1.2.6.
  \url{http://CRAN.R-project.org/package=fdrtool}

  R. Gentleman, V. Carey, W. Huber and F. Hahne.
  genefilter: methods for filtering genes from microarray experiments.
  R package version 1.28.2.
}

\keyword{package}

\seealso{
  \code{\link{gsri}}

  \code{\link{gsriFromFile}}
}

\examples{
## Simulate expression data for a gene set of:
## 100 genes, 40 samples (20 treatment, 20 control)
## and 30 regulated genes
expdata <- matrix(rnorm(4000,mean=0),nrow=100,ncol=40)
expdata[1:30,1:20] <- rnorm(600,mean=1)
data <- data.frame(expdata)
phenotype <- c(rep(0,20),rep(1,20))
geneSetName <- "Test Gene Set"

## Estimate the number of differentially expressed genes
res <- gsri(data, phenotype, geneSetName)


## Read expression data, phenotypes and gene sets from files
## Input files can be found in the vignette directory of this package
dataFileName <- system.file("extdata", "data.gct", package="GSRI")
phenotypeFileName <- system.file("extdata", "phenotype.cls", package="GSRI")
geneSetFileName <- system.file("extdata", "geneSets.gmt", package="GSRI")

res2 <- gsriFromFile(dataFileName, phenotypeFileName, geneSetFileName)
}
