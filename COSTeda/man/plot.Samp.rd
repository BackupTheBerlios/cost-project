\name{plot.Samp,Delta.list-method}
\alias{plot.Samp}
\alias{plot.Samp,Delta.list-method}
\docType{methods}
\title{Plot "Delta.length" object for a specified sample}
\description{
This method plots length distribution of one specified sample, compared with overall relative length distribution. Unlike \emph{plot.Delta.length}, it asks for an object of class \emph{Delta.list}.
}

\usage{
plot.Samp(x,SampNum,show.legend="right",\dots)
}

\arguments{
  \item{x}{A \emph{Delta.list} object created by \emph{Delta.cs} procedure.}
  \item{SampNum}{Character specifying sample Id.}
  \item{show.legend}{Display the legend (\code{""} means "no legend").}
  \item{...}{Further graphical arguments.}
}

\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{Delta.list}}, \code{\link{Delta.cs}}, \code{\link{Delta.length}}, \code{\link{plot.Delta.list}}, \code{\link{plot.Delta.length}}, \code{\link{plot.LD}}
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
plot.Samp(res,"194")
}
\keyword{methods}
