\name{ML.boxplot,csData-method}
\alias{ML.boxplot}
\alias{ML.boxplot,csData-method}
\docType{methods}
\title{Boxplot of individual biological parameters : Maturity by Length}
\description{
This method implements a boxplot of individual maturity-length data (maturity by length). It requires a \emph{csData} object with \emph{ca} table.
}

\usage{
ML.boxplot(object,\dots)
}

\arguments{
  \item{object}{A \emph{csData} object with \emph{ca} informations.}
  \item{...}{Further graphical arguments.}
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{WL.plot}}, \code{\link{ML.plot}}, \code{\link{SL.boxplot}}, \code{\link{L.plot.design}}
}

\examples{
data(sole)
ML.boxplot(sole.cs)
}
\keyword{methods}
