\name{plot.Delta.length-methods}
\docType{methods}
\alias{plot,Delta.length-method}
\alias{plot.Delta.length}
\alias{plot.Delta.length,Delta.length-method}
\title{Plot "Delta.length" object}
\description{
This method plots the output of outliers identification from \emph{plot.Delta.list} procedure as an object of class \emph{Delta.length}.
It consists in length distributions of each identified sample, compared with overall relative length distribution.
}

\usage{\S4method{plot}{Delta.length}(x,show.legend="right",\dots)}

\arguments{
  \item{x}{A \emph{Delta.length} object created by \emph{plot.Delta.list} procedure.}
  \item{show.legend}{Display the legend (\code{""} means "no legend").}
  \item{...}{Further graphical parameters.}
}

\section{Methods}{
\describe{
	\item{plot}{\code{signature(x="ANY",y="ANY")}: generic function: see 'plot'.}
	\item{plot}{\code{signature(Delta.length)}: plotting procedure of an object of class \emph{Delta.length}.}
}}


\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}

\author{Mathieu Merzereaud}
\seealso{\code{\link{Delta.length}}, \code{\link{Delta.list}}, \code{\link{Delta.cs}}, \code{\link{plot.Delta.list}}, \code{\link{plot.LD}}, \code{\link{plot.Samp}}, \code{\link{plot}}
}

\examples{
data(sole3.cs)
object <- sole3.cs
#only sea sampling data is kept
object@tr <- object@tr[object@tr$sampType=="S",]
object@hh <- object@hh[object@hh$sampType=="S",]
object@sl <- object@sl[object@sl$sampType=="S",]
object@hl <- object@hl[object@hl$sampType=="S",]

res <- Delta.cs(object,"SOL","LAN","quarter","month")
ident <- plot(res,selection=TRUE)
#plot(ident)
}
\keyword{dplot}
