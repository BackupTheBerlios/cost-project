\name{lenDisPlot,csData-method}
\alias{lenDisPlot}
\alias{lenDisPlot,csData-method}
\docType{methods}
\title{Plot of length distribution for a specified trip}
\description{
This method plots from a \emph{csData} object the length distribution of one or several samples in a specified trip, for given species and catch category.
}

\usage{
lenDisPlot(x,trpCode,species,fraction="LAN",staNum="all",\dots)
}

\arguments{
  \item{x}{A \emph{csData} object with \emph{hl} informations.}
  \item{trpCode}{Character specifying trip code.}
  \item{species}{Field specifying species (e.g \code{"Solea vulgaris"}).}                                                                                                                         
  \item{fraction}{Field specifying catch category (e.g \code{"LAN"}).}
  \item{staNum}{Character vector specifying sample(s) (or FOs). \code{"all"} results in displaying all samples of the specified trip, and \code{"allSum"} adds them up.}
  \item{...}{Further graphical arguments.}
}

\details{For more informations about arguments, see \emph{FishFrame/COST Exchange format specification}.}


\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{dltCls}}, \code{\link{dltId}}, \code{\link{dltCalc}}, \code{\link{dltPlot}}, \code{\link{plot.dltId}}, \code{\link{smpPlot}}
}

\examples{
data(sole)
lenDisPlot(sole.cs,"DIL1197","Solea vulgaris","DIS")
}
\keyword{methods}
