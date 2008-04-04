\name{L.plot.design,csData-method}
\alias{L.plot.design}
\alias{L.plot.design,csData-method}
\docType{methods}
\title{Plot.design of "csData" individual biological parameters in "ca"}
\description{
Version of \emph{plot.design} procedure from \pkg{graphics} package applied to \emph{csData} biological parameters.
}

\usage{
L.plot.design(object,\dots)
}

\arguments{
  \item{object}{A \emph{csData} object with \emph{ca} informations.}
  \item{...}{Further arguments.}
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{BioPar.plot}}, \code{\link{BioPar.boxplot}}
}

\examples{
data(sole)
L.plot.design(sole.cs)
}
\keyword{methods}
