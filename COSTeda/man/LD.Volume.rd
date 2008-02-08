\name{LD.Volume,csData-method}
\alias{LD.Volume}
\alias{LD.Volume,csData-method}
\docType{methods}
\title{Calculation upon sea-sampled catch at FO level, fishing day level, and trip level.}
\description{
This method creates an object of class \emph{LD.Vol} containing volume informations about sampled and raised catch,
for a given species in a specified catch category. It requires a \emph{csData} object built from \pkg{COSTcore} package.
}

\usage{
LD.Volume(object,fraction="all",species="all",TimeStrat=NULL,TechStrat=NULL,
         SpaceStrat=NULL,\dots)
}

\arguments{
  \item{object}{A \emph{csData} object with sea-sampling information (\emph{tr}, \emph{hh} and \emph{sl} required).}
  \item{fraction}{Field specifying catch category (to be chosen between \code{"LAN"}, \code{"DIS"} and \code{"all"}). See Details.}
  \item{species}{Field specifying species (e.g \code{"SOL"}). See Details.}
  \item{TimeStrat}{
  Specified time stratification (to be chosen between \code{"year"}, \code{"quarter"}, \code{"month"}
  and \code{NULL}). See Details.}
  \item{TechStrat}{
  Specified technical stratification (to be chosen between \code{"gear"}, \code{"foCatNat"}, \code{"foCatEu5"},
  \code{"foCatEu6"} and \code{NULL}). See Details.}
  \item{SpaceStrat}{
  Specified space stratification (to be chosen between \code{"area"}, \code{"rect"} and \code{NULL}). See Details.}
  \item{...}{Further arguments.}
}

\details{For more informations about arguments, see \emph{FishFrame/COST Exchange format specification}.}

\value{An object of class \emph{LD.Vol}.}

\author{Mathieu Merzereaud}

\seealso{\code{\link{LD.Vol}}, \code{\link{FD.plot}}, \code{\link{FD.boxplot}}, \code{\link{FO.plot}}, \code{\link{FO.boxplot}}
}

\examples{
data(sole3.cs)
object <- sole3.cs
#only sea sampling data is kept
object@tr <- object@tr[object@tr$sampType=="S",]
object@hh <- object@hh[object@hh$sampType=="S",]
object@sl <- object@sl[object@sl$sampType=="S",]
object@hl <- object@hl[object@hl$sampType=="S",]

x <- LD.Volume(object,fraction="LAN",species="SOL")
}

\keyword{methods}
