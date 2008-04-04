\name{plot.DeltaID-methods}
\docType{methods}
\alias{plot,DeltaID-method}
\alias{plot.DeltaID}
\alias{plot.DeltaID,DeltaID-method}
\title{Plot "DeltaID" object}
\description{
This method plots the output of outliers identification from \emph{plot.Delta} procedure as an object of class \emph{DeltaID}.
It consists in length distributions of each identified sample, compared with overall relative length distribution.
}

\usage{\S4method{plot}{DeltaID}(x,y=NULL,show.legend="right",\dots)}

\arguments{
  \item{x}{A \emph{DeltaID} object created by \emph{plot.Delta} procedure.}
  \item{y}{...}
  \item{show.legend}{Display the legend (\code{""} means "no legend").}
  \item{...}{Further graphical parameters.}
}

\section{Methods}{
\describe{
	\item{plot}{\code{signature(x="ANY",y="ANY")}: generic function: see 'plot'.}
	\item{plot}{\code{signature(DeltaID)}: plotting procedure of an object of class \emph{DeltaID}.}
}}


\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}

\author{Mathieu Merzereaud}
\seealso{\code{\link{DeltaID}}, \code{\link{DeltA-class}}, \code{\link{plot.Delta}}, \code{\link{Delta}}, \code{\link{plot.LD}}, \code{\link{plot.Samp}}, \code{\link{plot}}
}

\examples{
data(sole)
#obj <- plot(sole.cs,species="Solea vulgaris",tempStrata="quarter",techStrata="commCat",
#strat1="techStrata",strat2="tempStrata",strategy="cc",selection=TRUE,show.legend="right")
#plot(obj)
}
\keyword{dplot}
