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
  Estimates the number of differentially expressed genes in gene sets.
}

\usage{
gsri(exprs, phenotype, geneSet, names=NULL, weight=NULL, nBoot=100,
test=rowt, testArgs=NULL, alpha=0.05, grenander=TRUE)
}

\arguments{
  \item{exprs}{A matrix or an expression set of size n x m, containing
    the expression values for the genes, with rows representing the
    genes and columns representing the arrays.}
  \item{phenotype}{A factor of length m representing the phenotypes of
    the samples.}
  \item{geneSet}{missing, GeneSet, GeneSetClass}
  \item{names}{Character vector with the names of the gene set(s).}
  \item{weight}{Weights}
  \item{nBoot}{An integer specifying the number of bootstrap samples to
    be drawn (default: 100).}
  \item{test}{A function specifying the statistical test used to assess
    the effect between the groups in terms of p-values. The groups are
    given by the 'phenotype' argument. Here, a t-test (rowt) and F-test
    (rowF) are implemented (default: rowt). A user-defined function can
    also be given, in order to specify a custom test statistics. For
    details see the 'Notes' section.}
  \item{testArgs}{Optional arguments passed to the function 'test' if
    specified (default: NULL)}
  \item{alpha}{A numeric specifying the significance level for
    bootstrapping. The resulting GSRI is the lower bound of the
    (1-\eqn{alpha}{alpha})*100\% confidence interval.}
  \item{grenander}{A logical indicating whether the modified Grenander
    estimator for the cumulative density should be used, instead of a
    classical ECDF (default: TRUE).}
}

\section{Methods}{
  \describe{
    \item{
      \code{signature(exprs="matrix", phenotype="factor", geneSet="missing")}
    }{bla}

    \item{
      \code{signature(exprs="ExpressionSet", phenotype="factor", geneSet="missing")}
    }{bla}

    \item{
      \code{signature(exprs="matrix", phenotype="factor", geneSet="GeneSet")}
    }{bla}
    
    \item{
      \code{signature(exprs="ExpressionSet", phenotype="factor", geneSet="GeneSet")}
    }{bla}
    
    \item{
      \code{signature(exprs="matrix", phenotype="factor", geneSet="GeneSetCollection")}
    }{bla}
    
    \item{
      \code{signature(exprs="ExpressionSet", phenotype="factor", geneSet="GeneSetCollection")}
    }{bla}
    }
}

\details{
  This function calculates the number of differentially expressed genes
  for a single gene set, with data and gene set taken from the
  workspace.
  
  From bootstrapping the group samples the (1-alpha)*100\% quantile of
  the distribution of the estimated number of differentially expressed
  genes is obtained. The Gene Set Regulation Index (GSRI) is defined as
  the 5\% quantile and indicates, that with a probability of 95\% more
  than GSRI genes in the gene set are differentially expressed.

  This index can be employed to test the hypothesis whether at least one
  gene in a set is regulated and to compare and rank the regulation of
  different gene sets.
}

\value{
  An object of class 'Gsri' with the slots:

  \item{result}{Result of the estimation.}
  \item{pval}{List with the p-values for each gene set.}
  \item{cdf}{List with the cumulative denstiy for each gene set.}
  \item{parms}{Parameter values used in the analysis.}
}

\author{
  Kilian Bartholome, Julian Gehring
}

\note{
  Usage of the Grenander estimate is based on the assumption about the
  concave shape of the cumulative density distribution. It reduces the
  variance, i.e.  makes the approach more stable especially for small
  gene sets. On the downside the number of significant genes is
  overestimated for few regulated genes. A conservative solution of this
  trade-off would be to rank the gene-sets with and without Grenander
  estimate and to choose the lowest rank for each gene-set.  Please note
  that the Grenander estimate does not allow duplicates in the p-values.
  If this occurs in a data set, an error message will be printed and the
  analysis should be repeated without the Grenander estimate.

  With the t-test and the F-test, two widely used statistical tests are
  available in this package. To allow fast computation this package uses
  the implementations from the \pkg{genefilter} package.

  It is also possible to apply user-defined tests with this method. In
  this case the function has to be called by function(data, phenotype,
  testArgs). 'data' and 'phenotype' are of class 'matrix' and 'factor',
  respectively. 'testArgs' will only be passed to the function if it is
  defined. In general all methods that compute p-values are suitable. The
  function has to return a vector with one p-value per gene. For details
  on how to use your own test functions please refer to the vignette of
  this package.
}

\seealso{
  \code{\link{GSRI}}
}

\examples{
\dontrun{
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
}
}
\keyword{distribution}
\keyword{htest}
