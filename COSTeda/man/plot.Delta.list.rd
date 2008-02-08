\name{plot.Delta.list-methods}
\docType{methods}
\alias{plot,Delta.list-method}
\alias{plot.Delta.list}
\alias{plot.Delta.list,Delta.list-method}
\title{Plot "Delta.list" object}
\description{
This method plots the output of \emph{Delta.cs} procedure as an object of class \emph{Delta.list}. It allows outliers identification to create a new object of class \emph{Delta.length}.
}

\usage{\S4method{plot}{Delta.list}(x,selection=FALSE,show.legend="right",\dots)}

\arguments{
  \item{x}{A \emph{Delta.list} object created by \emph{Delta.cs} procedure.}
  \item{selection}{If \code{TRUE}, outliers identification is made, and a \emph{Delta.length} object is returned. 
  Displayed numbers during identification process are numbers of corresponding \emph{hh} dataset lines.}
  \item{show.legend}{Display the legend (\code{""} means "no legend").}
  \item{...}{Further graphical parameters.}
}

\section{Methods}{
\describe{
	\item{plot}{\code{signature(x="ANY",y="ANY")}: generic function: see 'plot'.}
	\item{plot}{\code{signature(Delta.list)}: plotting procedure of an object of class \emph{Delta.list}.}
}}


\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}

\author{Mathieu Merzereaud}
\seealso{\code{\link{Delta.list}}, \code{\link{Delta.cs}}, \code{\link{Delta.length}}, \code{\link{plot.Delta.length}}, \code{\link{plot.LD}}, \code{\link{plot.Samp}}, \code{\link{plot}}
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
#select FOs (points) and press right-click to end
}

\keyword{dplot}
