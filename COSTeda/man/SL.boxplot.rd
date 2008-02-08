\name{SL.boxplot,csData-method}
\alias{SL.boxplot}
\alias{SL.boxplot,csData-method}
\docType{methods}
\title{Boxplot of individual biological parameters : Sex by Length}
\description{
This method implements a boxplot of individual sex-length data (sex by length). It requires a \emph{csData} object with \emph{ca} table.
}

\usage{
SL.boxplot(object,\dots)
}

\arguments{
  \item{object}{A \emph{csData} object with \emph{ca} informations.}
  \item{...}{Further graphical arguments.}
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{WL.plot}}, \code{\link{ML.plot}}, \code{\link{ML.boxplot}}, \code{\link{L.plot.design}}
}

\examples{
data(sole3.cs)
SL.boxplot(sole3.cs)
}
\keyword{methods}
