\name{GSRI-package}

\docType{package}

\alias{GSRI}
\alias{GSRI-package}

\title{Gene Set Regulation Index (GSRI) package}

\description{
  The \pkg{GSRI} package estimates the number of differentially
  expressed genes in gene sets, utilizing the concept of the Gene Set
  Regulation Index (GSRI).
}

\details{
  The GSRI approach estimates the number of differentially expressed
  genes in gene sets. It does not require a cut-off value for the
  distinction between regulated and unregulated genes. The approach is
  based on the fact that p-values obtained from a statistical test are
  uniformly distributed under the null hypothesis and shifted towards
  zero in case of the alternative hypothesis.

  The Gene Set Regulation Index (GSRI) is the
  \eqn{alpha*100\%}{\alpha*100\%}\% quantile of the distribution of the
  estimated number of differentially expressed genes, obtained from
  bootstrapping the arrays within the groups. It indicates that with a
  probability of \eqn{(1-alpha)*100\%}{(1-\alpha)*100\%} more than
  \eqn{eta}{\eta} genes in the gene set are differentially
  expressed. This index can be employed to rank and compare the
  regulation of gene sets as well as to test the hypothesis whether at
  least one gene in the gene set is regulated.
  
  For details of the method, an application to experimental data, and a
  comparison with related approaches, please see Bartholomé et al.,
  2009.

  The package is published under the GPL-3 license.
}

\author{
  Julian Gehring, Kilian Bartholomé, Clemens Kreutz, Jens Timmer

  Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}

\references{
  Kilian Bartholomé, Clemens Kreutz, and Jens Timmer:
  Estimation of gene induction enables a relevance-based ranking of gene sets,
  Journal of Computational Biology: A Journal of Computational Molecular Cell 
  Biology 16, no. 7 (July 2009): 959-967.
  \url{http://www.liebertonline.com/doi/abs/10.1089/cmb.2008.0226}

  The \pkg{GSRI} package uses the functionality of the following
  packages:

  Julian Gehring, Clemens Kreutz, Jens Timmer:
  les: Identifying Loci of Enhanced Significance in Tiling Microarray Data
  \url{http://bioconductor.org/help/bioc-views/release/bioc/html/les.html}

  Korbinian Strimmer:
  fdrtool: Estimation and Control of (Local) False Discovery Rates.
  \url{http://CRAN.R-project.org/package=fdrtool}

  Robert Gentleman, Vincent J. Carey, Wolfgang Huber, Florian Hahne:
  genefilter: methods for filtering genes from microarray experiments.
  \url{http://bioconductor.org/help/bioc-views/release/bioc/html/genefilter.html}
}

\seealso{
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

\dontrun{
## Read expression data, groups and gene sets from files
## Input files can be found in the vignette directory of this package
dataFileName <- system.file("extdata", "data.gct", package="GSRI")
groupsFileName <- system.file("extdata", "phenotype.cls", package="GSRI")
geneSetFileName <- system.file("extdata", "geneSets.gmt", package="GSRI")

res2 <- gsriFromFile(dataFileName, groupsFileName, geneSetFileName)
}
}

\keyword{package}
\keyword{htest}
\keyword{distribution}