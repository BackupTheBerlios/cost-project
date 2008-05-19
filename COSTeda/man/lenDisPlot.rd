\name{lenDisPlot}
\alias{lenDisPlot}
\alias{lenDisPlot,csData-method}
\alias{lenDisPlot,csDataVal-method}
\docType{methods}
\title{Plot of length distribution for a specified trip}
\description{
This method plots from a \emph{csData/csDataVal} object the length distribution of one or several samples in a specified trip, for given species and catch category.
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

\author{Mathieu Merzereaud}

\seealso{\code{\link{edaResult}}, \code{\link{deltCalc}}, \code{\link{plot.edaResult}}
}

\examples{
data(sole)
lenDisPlot(sole.cs,"DIL1197","Solea vulgaris","DIS")
}
\keyword{methods}
