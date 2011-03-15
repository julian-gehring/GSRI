\name{plot-methods}
\docType{methods}

\alias{plot}
\alias{plot,Gsri,ANY-method}
\alias{plot,Gsri,integer-method}

\title{
  Visualize GSRI results
}

\description{
  Plot the empirical cumulative density along with the estimated degree
  of regulation of the p-value distribution for a gene set.
}

\usage{
plot(x, y, ...)
}

\arguments{
  \item{x}{An object of class 'Gsri' containing the results to plot.}
  \item{y}{A single integer or character string specifying the results
    of which gene set to plot. This has to be given if 'x' contains the
    results of several gene sets. An integer is interpreted as the index
    of the gene set (i.e. the row number), while a character is matched
    against the names of the gene sets.}
  \item{...}{Optional arguments used in order to customize the plot. See
  the 'details' section.}
}

\details{
  Plotting the p-value distribution and the estimation of the
  regularized compontent for a gene set allows to insepct the results
  in detail. The plot illustrates the empirical cumulative density
  function of the p-values obtained from testing for a differential
  effect between the groups. Additionally, the fit of the uniformly
  distributed component along with the estimated fraction of regularized
  genes and the GSRI is shown.

  The 'plot' method uses a special system in order to customize the
  graphical elements of the figure. It allows to refer to the different
  components with the name of the additional input argument; its value
  is a list containing named graphical parameters for the underlying
  plot function. The following list describes the possible names and
  their contribution.

  \describe{
    \item{'plot'}{Arguments for the axis and the labeling, passed to the
      'plot' function.}
    \item{'fit'}{Arguments for the fit of the linear component of the
      ECDF, corresponding to the part without differential effect,
      passed to the 'lines' function.}
    \item{'ecdf'}{Arguments for the ECDF of the p-values, passed to the
      'lines' function.}
    \item{'reg'}{Arguments for the horizontal line indicating the fraction
      of regulation, passed to the 'lines' function.}
    \item{'regText'}{Arguments for the text label of 'reg', passed to
      the 'text' function.}
    \item{'gsri'}{Arguments for the horizontal line indicating the
      GSRI, passed to the 'lines' function.}
    \item{'gsriText'}{Arguments for the text label of 'gsri', passed to
      the 'text' function.}
  }

  Thus, changing for example the limit of the y-axis, the plot type and
  color of the ECDF, and the label of the x-axis, you can use:

  \code{plot(x, plot=list(ylim=c(0, 0.8), xlab=expression(p)),
  ecdf=list(type="s", col="darkgreen"))}

  For more details, please see the 'examples' section.
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
\dontrun{
plot(x)

plot(x, plot=list(main="Example plot"), ecdf=list(pch=21),
fit=list(lty=2, lwd=0.5, col="black"), gsri=list(col="lightblue"))

plot(x2, 2)
plot(x2, "gs2")
}
}

\keyword{utilities}
\keyword{methods}
\keyword{hplot}
