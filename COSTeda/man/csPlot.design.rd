\name{csPlot.design,csDataVal-method}
\alias{csPlot.design}
\alias{csPlot.design,csData-method}
\alias{csPlot.design,csDataVal-method}
\docType{methods}
\title{Plot.design of individual biological parameters in "ca" table from a "csData" or "csDataVal" object}
\description{
Version of \emph{plot.design} procedure from \pkg{graphics} package applied to \emph{csData} or \emph{csDataVal} biological parameters.
}

\usage{
csPlot.design(object,species="all",\dots)
}

\arguments{
  \item{object}{A \emph{csData} or \emph{csDataVal} object with \emph{ca} informations.}
  \item{species}{A character (vector) specifying species (default value is "all" for all species) .}
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
