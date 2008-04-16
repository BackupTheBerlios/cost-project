\name{csPlot.design,csDataVal-method}
\alias{csPlot.design}
\alias{csPlot.design,csDataVal-method}
\docType{methods}
\title{Plot.design of individual biological parameters in "ca" table from "csDataVal" object}
\description{
Version of \emph{plot.design} procedure from \pkg{graphics} package applied to \emph{csDataVal} biological parameters.
}

\usage{
csPlot.design(object,\dots)
}

\arguments{
  \item{object}{A \emph{csDataVal} object with \emph{ca} informations.}
  \item{...}{Further arguments.}
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{bioPar.plot}}, \code{\link{bioPar.boxplot}}
}

\examples{
data(sole)
sole.cs.val <- csDataVal(sole.cs)
csPlot.design(sole.cs.val)
}
\keyword{methods}
