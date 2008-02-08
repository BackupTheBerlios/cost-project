\name{plot.LD,csData-method}
\alias{plot.LD}
\alias{plot.LD,csData-method}
\docType{methods}
\title{Plot of length distribution for a specified trip}
\description{
This method plots from a \emph{csData} object the length distribution of one or several samples in a specified trip, for given species and catch category.
}

\usage{
plot.LD(x,species,fraction="LAN",trpCode,staNum="all",\dots)
}

\arguments{
  \item{x}{A \emph{csData} object with \emph{hl} informations.}
  \item{species}{Field specifying species (e.g \code{"SOL"}). See Details.}
  \item{fraction}{Field specifying catch category (e.g \code{"LAN"}). See Details.}
  \item{trpCode}{Character specifying trip code. See Details.}
  \item{staNum}{Character vector specifying sample(s) (or FOs). See Details.}
  \item{...}{Further graphical arguments.}
}

\details{For more informations about arguments, see \emph{FishFrame/COST Exchange format specification}.}


\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{Delta.list}}, \code{\link{Delta.cs}}, \code{\link{Delta.length}}, \code{\link{plot.Delta.list}}, \code{\link{plot.Delta.length}}, \code{\link{plot.Samp}}
}

\examples{
data(sole3.cs)
#length distribution by trip
plot.LD(sole3.cs,"SOL","LAN","LIM1")      #sampType="S"
plot.LD(sole3.cs,"SOL","LAN","131395")    #sampType="M"

#for all the trip
temp <- sole3.cs
temp@hl$staNum <- "999"
plot.LD(temp,"SOL","LAN","LIM1")
}
\keyword{methods}
