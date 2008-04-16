\name{plot.dltId-methods}
\docType{methods}
\alias{plot,dltId-method}
\alias{plot.dltId}
\alias{plot.dltId,dltId-method}
\title{Plot "dltId" object}
\description{
This method plots the output of outliers identification from \emph{dltPlot} procedure as an object of class \emph{dltId}.
It consists in length distributions of each identified sample, compared with overall relative length distribution.
}

\usage{\S4method{plot}{dltId}(x,y=NULL,show.legend="right",\dots)}

\arguments{
  \item{x}{A \emph{dltId} object created by \emph{dltPlot} procedure.}
  \item{y}{...}
  \item{show.legend}{Display the legend (\code{""} means "no legend").}
  \item{...}{Further graphical parameters.}
}

\section{Methods}{
\describe{
	\item{plot}{\code{signature(x="ANY",y="ANY")}: generic function: see 'plot'.}
	\item{plot}{\code{signature(dltId)}: plotting procedure of an object of class \emph{dltId}.}
}}


\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}

\author{Mathieu Merzereaud}
\seealso{\code{\link{dltCls}}, \code{\link{dltId}}, \code{\link{dltCalc}}, \code{\link{dltPlot}}, \code{\link{smpPlot}}, \code{\link{lenDisPlot}}
}

\examples{
data(sole)
#obj <- dltPlot(sole.cs,species="Solea vulgaris",timeStrata="quarter",techStrata="commCat",
#strat1="techStrata",strat2="timeStrata",strategy="cc",selection=TRUE,show.legend="right")
#plot(obj)
}
\keyword{dplot}
