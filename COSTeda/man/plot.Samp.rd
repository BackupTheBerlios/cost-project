\name{plot.Samp,DeltaID-method}
\alias{plot.Samp}
\alias{plot.Samp,DeltaID-method}
\docType{methods}
\title{Plot "DeltaID" object for a specified sample}
\description{
This method plots length distribution of one specified sample, compared with overall relative length distribution.
}         

\usage{
plot.Samp(x,SampNum,show.legend="right",\dots)
}

\arguments{
  \item{x}{A \emph{DeltaID} object created by \emph{plot.Delta} procedure.}
  \item{SampNum}{Character specifying sample Id.}
  \item{show.legend}{Display the legend (\code{""} means "no legend").}
  \item{...}{Further graphical arguments.}
}

\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{DeltaID}}, \code{\link{DeltA-class}}, \code{\link{plot.Delta}}, \code{\link{plot.DeltaID}}, \code{\link{plot.LD}}, \code{\link{Delta}}
}

\examples{
data(sole)
#obj2 <- plot(sole.cs,species="Solea vulgaris",tempStrata="quarter",techStrata="commCat",strat1="techStrata",strat2="tempStrata",selection=TRUE,show.legend="right")
##"xxx" is the index of an identified sample
#plot.Samp(obj2,"xxx") 
}
\keyword{methods}
