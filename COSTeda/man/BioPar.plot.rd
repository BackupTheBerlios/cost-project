\name{BioPar.plot,csData-method}
\alias{BioPar.plot}
\alias{BioPar.plot,csData-method}
\docType{methods}
\title{Plots of individual biological parameters}
\description{
This method implements a scatter plot of individual biological data. It requires a \emph{csData} object with \emph{ca} table.
}

\usage{
BioPar.plot(object,type="WL",selection=FALSE,\dots)}

\arguments{
  \item{object}{A \emph{csData} object with \emph{ca} informations.}
  \item{type}{To be chosen between \code{"WL"} (weight~length), \code{"ML"} (maturity~length) or \code{"SL"} (sex~length).}
  \item{selection}{If \code{TRUE}, outlier(s) identification can be done and corresponding data is returned.}
  \item{...}{Further graphical arguments.}
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{BioPar.boxplot}}, \code{\link{L.plot.design}}
}

\examples{
data(sole)
BioPar.plot(sole.cs)
}
\keyword{methods}
