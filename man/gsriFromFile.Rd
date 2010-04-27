\name{gsriFromFile}

\alias{gsriFromFile}

\title{
  Gene Set Regulation Index (GSRI)
}

\description{
  Estimates the number of differentially expressed genes for a list of
  gene sets, with data and gene sets read from files.
}

\usage{
gsriFromFile(dataFileName, phenotypeFileName, geneSetFileName,
             useGrenander = FALSE, plotResults = FALSE, writeResults =
             FALSE, writeSummary = FALSE, minGeneSetSize = 10,
             nBootstraps = 100, test = "ttest", testArgs = NULL, alpha =
             0.05, verbose = TRUE, prec = 4)
}

\arguments{
  \item{dataFileName}{Data file of type '.gct' containing the expression
    data. For details on the format see the Notes section below.}
  \item{phenotypeFileName}{Phenotype file of type '.cls' specifying the
    phenotypes of the samples. For details on the format see the Notes
    section below.}
  \item{geneSetFileName}{Gene set file of type '.gmt' containing a list
    of gene sets. For details on the format see the Notes section
    below.}
  \item{useGrenander}{Logical indicating whether the grenander estimate
    from the \pkg{fdrtool} package should be used additionally (default:
    FALSE). For details see the Notes section below.}
  \item{plotResults}{Logical indicating whether the results should be
    plotted (default: FALSE). The plot shows the cumulative density
    function of calculated p-values, the fit of the uniform
    distribution, the estimated number of regulated genes and the
    estimated GSRI for each gene set.}
  \item{writeResults}{Logical indicating whether a list of the estimated
    regulated genes for each gene set should be written to a text file
    in the working directory (default: FALSE). The file will be named
    'GeneSet_#geneSetName#_data.txt', with #geneSetName# taken from the
    'geneSetName' argument.}
  \item{writeSummary}{Logical indicating whether a summary of results
    for all gene sets should be written into a single file (default:
    FALSE). The file will be named '#dataFileName#_results.txt', with
    #dataFileName# taken from the file name of the 'dataFileName'
    argument.}
  \item{minGeneSetSize}{Integer specifying the minimal gene set size
    (default: 10). Only gene sets having at least this size will be used
    for estimation.}
  \item{nBootstraps}{Number of bootstrap samples to be drawn (default:
    100)}
  \item{test}{Character string or function name to specify the
    statistical test to calculate p-values for effect between
    groups. Groups are defined by the 'phenotype' argument. In this
    package both a t-test (default: "ttest") and an F-test ("ftest")
    between groups are implemented and can be chosen with the according
    character string. An user-defined function can also be passed as an
    argument in order to specify own test statistics. For details see
    the Notes section below.}
  \item{testArgs}{Optional arguments passed to the function 'test' if
    specified (default: NULL)}
  \item{alpha}{Significance level for bootstrap (default: 0.05). The
    resulting GSRI will be the (1-alpha)*100\% confidence interval. If a
    vector of values is given, GSRI will be calculated for all values
    and passed to the output argument. Plots and file outputs will only
    contain GSRI for the first value in the vector to simplifyoutput.}
  \item{verbose}{Logical indicating whether information about data
    loading and gene set computation should be printed on screen
    (default: TRUE).}
  \item{prec}{Integer specifing the number of signifiant digits to write
    to files. For details see 'signif'.}
}

\details{
  This function calculates the number of differentially expressed genes
  for a list of gene sets, with data and gene sets read from files.

  From bootstrapping the group samples the (1-alpha)*100\% quantile of the
  distribution of the estimated number of differentially expressed genes
  is obtained. The Gene Set Regulation Index (GSRI) is defined as the 5\%
  quantile and indicates, that with a probability of 95\% more than GSRI
  genes in the gene set are differentially expressed.

  This index can be employed to test the hypothesis whether at least one
  gene in a set is regulated and to compare and rank the regulation of
  different gene sets.
}

\value{
  A list with components (one entry per gene set): 
   
  \item{geneSet}{Names of gene sets}
  \item{percRegGenes}{Estimated percentage of differentially expressed
    genes}
  \item{percRegGenesSd}{Estimated standard deviation of the percentage
    of differentially expressed genes}
  \item{numRegGenes}{Estimated number of differentially expressed
    genes}
  \item{numRegGenesSd}{Estimated standard deviation of the number of
    differentially expressed genes}
  \item{gsri}{Gene Set Regulation Index}
  \item{nGenes}{Number of genes in gene sets}
  \item{geneSetGenes}{Names of genes in gene sets.}
}

\author{
  Kilian Bartholomé, Julian Gehring
}

\note{
  The input files should have the format typical for those file
  types. Details on the specific formats can be found at
  \url{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}.
  For more details and example files please refer to the vignette of
  this package.

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
  defined. In general all methods that compute p-values are
  suitable. The function has to return a vector with one p-value per
  gene. For details on how to use your own test functions please refer
  to the vignette of this package.
}

\seealso{
\code{\link{gsri}}

\code{\link{GSRI}}
}

\examples{
## Read expression data, phenotypes and gene sets from files
## Input files can be found in the vignette directory of this package
dataFileName <- system.file("extdata", "data.gct", package="GSRI")
phenotypeFileName <- system.file("extdata", "phenotype.cls", package="GSRI")
geneSetFileName <- system.file("extdata", "geneSets.gmt", package="GSRI")

res <- gsriFromFile(dataFileName, phenotypeFileName, geneSetFileName)
}

\keyword{distribution}
