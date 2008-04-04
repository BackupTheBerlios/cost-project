\name{LD.Volume,csData-method}
\alias{LD.Volume}
\alias{LD.Volume,csData-method}
\docType{methods}
\title{Calculation upon sea-sampled catch at FO level, fishing day level, and trip level.}
\description{
This method creates an object of class \emph{LD.Vol} containing volume informations about sampled and raised catch,
for a given species in a specified catch category. It requires a \emph{csData} object built from \pkg{COSTcore} package.
Only sea sampling data is computed.
}

\usage{
LD.Volume(object,species,fraction="LAN",TimeStrat=NULL,TechStrat=NULL,
         SpaceStrat=NULL,\dots)
}

\arguments{
  \item{object}{A \emph{csData} object with sea-sampling information (\emph{tr}, \emph{hh} and \emph{sl} required).}
  \item{species}{Field specifying species (e.g \code{"Solea vulgaris"}). See Details.}
  \item{fraction}{Field specifying catch category (to be chosen between \code{"LAN"} and \code{"DIS"}.). See Details.}
  \item{TimeStrat}{
  Specified time stratification (to be chosen between \code{"year"}, \code{"quarter"}, \code{"month"}, \code{NULL},...). See Details.}
  \item{TechStrat}{
  Specified technical stratification (to be chosen between \code{"gear"}, \code{"foCatNat"}, \code{"foCatEu5"}, \code{NULL},...). See Details.}
  \item{SpaceStrat}{
  Specified space stratification (to be chosen between \code{"area"}, \code{"rect"}, \code{NULL},...). See Details.}
  \item{...}{Further arguments.}
}

\details{For more informations about arguments, see \emph{FishFrame/COST Exchange format specification}.}

\value{An object of class \emph{LD.Vol}.}

\author{Mathieu Merzereaud}

\seealso{\code{\link{LD.Vol}}, \code{\link{FD.plot}}, \code{\link{FD.boxplot}}, \code{\link{FO.plot}}, \code{\link{FO.boxplot}}
}

\examples{
data(sole)
object <- sole.cs

x <- LD.Volume(object,fraction="LAN",species="Solea vulgaris")
}

\keyword{methods}
