\name{smpPlot,dltId-method}
\alias{smpPlot}
\alias{smpPlot,dltId-method}
\docType{methods}
\title{Plot "dltId" object for a specified sample}
\description{
This method plots length distribution of one identified sample, compared with overall relative length distribution.
}         

\usage{
smpPlot(x,smpNum,show.legend="right",\dots)
}

\arguments{
  \item{x}{A \emph{dltId} object created by \emph{dltPlot} procedure.}
  \item{smpNum}{Character specifying sample Id as displayed during outliers identification process.}
  \item{show.legend}{Display the legend (\code{""} means "no legend").}
  \item{...}{Further graphical arguments.}
}

\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{dltCls}}, \code{\link{dltId}}, \code{\link{dltCalc}}, \code{\link{dltPlot}}, \code{\link{plot.dltId}}, \code{\link{lenDisPlot}}
}

\examples{
data(sole)
#obj <- dltPlot(sole.cs,species="Solea vulgaris",timeStrata="quarter",techStrata="commCat",
#strat1="techStrata",strat2="timeStrata",strategy="cc",selection=TRUE,show.legend="right")
##'xxx' is the index of an identified sample
#smpPlot(obj,"xxx") 
}
\keyword{methods}
