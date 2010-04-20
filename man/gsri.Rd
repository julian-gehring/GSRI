\name{gsri}
\alias{gsri}
\title{
Gene Set Regulation Index (GSRI)
}
\description{
Estimates the number of differentially expressed genes for a single gene set.
}
\usage{
gsri(data, phenotype, geneSetName, 
     useGrenander = FALSE, plotResults = TRUE, writeResults = FALSE,
     nBootstraps = 100, test = "ttest", testArgs = NULL, alpha = 0.05)
}

\arguments{
  \item{data}{A data frame or matrix of size n x m containing the gene expression dataset with rows representing genes and columns samples.}
  \item{phenotype}{Vector of size m containing the phenotypes of the samples. If not already a factorial variable it will be internally converted to a factor.}
  \item{geneSetName}{Character string with the name of the gene set. Used for identification of the data set, has no influence on the calculation.}
  \item{useGrenander}{Logical indicating whether the grenander estimate from the \pkg{fdrtool} package should be used additionally (default: FALSE). For details see the Notes section below.}
  \item{plotResults}{Logical indicating whether the results should be plotted (default: FALSE). The plot shows the cumulative density function of calculated p-values, the fit of the uniform distribution, the estimated number of regulated genes and the estimated GSRI for each gene set.}
  \item{writeResults}{Logical indicating whether a list of the estimated regulated genes for each gene set should be written to a text file in the working directory (default: FALSE). The file will be named 'GeneSet_#geneSetName#_data.txt', with #geneSetName# taken from the 'geneSetName' argument.}
  \item{nBootstraps}{Number of bootstrap samples to be drawn (default: 100)}
  \item{test}{Character string or function name to specify the statistical test to calculate p-values for effect between groups. Groups are defined by the 'phenotype' argument. In this package both a t-test (default: "ttest") and an F-test ("ftest") between groups are implemented and can be chosen with the according character string. A user-defined function can also be passed as an argument in order to specify own test statistics. For details see the Notes section below.}
  \item{testArgs}{Optional arguments passed to the function 'test' if specified (default: NULL)}
  \item{alpha}{Significance level for bootstrap (default: 0.05). The resulting GSRI will be the (1-alpha)*100\% confidence interval. If a vector of values is given, GSRI will be calculated for all values and passed to the output argument. Plots and file outputs will only contain GSRI for the first value in the vector to simplify output.}
}
\details{
This function calculates the number of differentially expressed genes for a 
single gene set, with data and gene set taken from the workspace.

From bootstrapping the group samples the (1-alpha)*100\% quantile of the 
distribution of the estimated number of differentially expressed genes is 
obtained. The Gene Set Regulation Index (GSRI) is defined as the 5\% quantile and
indicates, that with a probability of 95\% more than GSRI genes in the gene set 
are differentially expressed.

This index can be employed to test the hypothesis whether at least one gene 
in a set is regulated and to compare and rank the regulation of different gene sets.
}
\value{
   A list with components:
   
  	\item{geneSet}{Name of gene set}
  	\item{percRegGenes}{Estimated percentage of differentially expressed genes}
  	\item{percRegGenesSd}{Estimated standard deviation of the percentage of differentially expressed genes}
  	\item{numRegGenes}{Estimated number of differentially expressed genes}
  	\item{numRegGenesSd}{Estimated standard deviation of the number of differentially expressed genes}
  	\item{gsri}{Gene Set Regulation Index}
  	\item{nGenes}{Number of genes in gene set}
}
\author{
Kilian Bartholomé, Julian Gehring
}
\note{
Usage of the Grenander estimate is based on the assumption about the concave 
shape of the cumulative density distribution. It reduces the variance, i.e. 
makes the approach more stable especially for small gene sets. On the downside 
the number of significant genes is overestimated for few regulated genes. A 
conservative solution of this trade-off would be to rank the gene-sets with 
and without Grenander estimate and to choose the lowest rank for each gene-set. 
Please note that the Grenander estimate does not allow duplicates in the p-values. 
If this occurs in a data set, an error message will be printed and the analysis 
should be repeated without the Grenander estimate.

With the t-test and the F-test, two widely used statistical tests are available 
in this package. To allow fast computation this package uses the 
implementations from the \pkg{genefilter} package.

It is also  possible to apply user-defined tests with this method. In this case the function 
has to be called by function(data, phenotype, testArgs). 'data' and 
'phenotype' are of class 'matrix' and 'factor', respectively. 'testArgs' will only
be passed to the function if it is defined. In general all methods  that compute 
p-values are suitable. The function has to return a vector with one 
p-value per gene. For details on how to use your own test functions please 
refer to the vignette of this package.
}

\seealso{
\code{\link{gsriFromFile}}

\code{\link{GSRI}}
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
}

\keyword{distribution}
