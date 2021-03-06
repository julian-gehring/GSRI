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
  genes in gene sets. It is independent of the underlying statistical
  test used for assessing the differential effect of genes and does not
  require any cut-off values for the distinction between regulated and
  unregulated genes. The approach is based on the fact that p-values
  obtained from a statistical test are uniformly distributed under the
  null hypothesis and are shifted towards zero in case of the
  alternative hypothesis.

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
  
  For details of the method, an application to experimental data, and a
  comparison with related approaches, see Bartholome et al., 2009.

  The package is published under the GPL-3 license.
}

\author{
  Julian Gehring, Kilian Bartholome, Clemens Kreutz, Jens Timmer

  Maintainer: Julian Gehring <julian.gehring@fdm.uni-freiburg.de>
}

\references{
  Kilian Bartholome, Clemens Kreutz, and Jens Timmer:
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
rownames(exprs) <- paste("g", 1:nrow(exprs), sep="")
groups <- factor(rep(1:2, each=10))

## Estimate the number of differentially expressed genes
res <- gsri(exprs, groups)
res

## Perform the analysis for different gene set
library(GSEABase)
gs1 <- GeneSet(paste("g", 25:40, sep=""), setName="set1")
gs2 <- GeneSet(paste("g", seq(1, nrow(exprs), by=5), sep=""), setName="set2")
gsc <- GeneSetCollection(gs1, gs2)

res2 <- gsri(exprs, groups, gs1)
res3 <- gsri(exprs, groups, gsc, verbose=TRUE)

summary(res2)
}

\keyword{package}
\keyword{htest}
\keyword{distribution}