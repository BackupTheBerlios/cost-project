\name{ML.plot,csData-method}
\alias{ML.plot}
\alias{ML.plot,csData-method}
\docType{methods}
\title{Plots of individual biological parameters : Maturity-Length}
\description{
This method implements a scatter plot of individual maturity-length data. It requires a \emph{csData} object with \emph{ca} table.
}

\usage{
ML.plot(object,\dots)
}

\arguments{
  \item{object}{A \emph{csData} object with \emph{ca} informations.}
  \item{...}{Further graphical arguments.}
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{WL.plot}}, \code{\link{ML.boxplot}}, \code{\link{SL.boxplot}}, \code{\link{L.plot.design}}
}

\examples{
data(sole)
ML.plot(sole.cs)
}
\keyword{methods}
