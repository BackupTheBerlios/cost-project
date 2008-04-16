\name{landisVolume,csData-method}
\alias{landisVolume}
\alias{landisVolume,csData-method}
\docType{methods}
\title{Calculation upon sea-sampled catch at FO level, fishing day level, and trip level.}
\description{
This method creates an object of class \emph{landisVol} containing volume informations about sampled and raised catch,
for a given species in a specified catch category. It requires a \emph{csData} object built from \pkg{COSTcore} package.
Only sea sampling data is computed.
}

\usage{
landisVolume(object,species,fraction="LAN",timeStrata=NULL,spaceStrata=NULL,
         techStrata=NULL,\dots)
}

\arguments{
  \item{object}{A \emph{csData} object with sea-sampling information (\emph{tr}, \emph{hh} and \emph{sl} required).}
  \item{species}{Field specifying species (e.g \code{"Solea vulgaris"}).}
  \item{fraction}{Field specifying catch category (to be chosen between \code{"LAN"} and \code{"DIS"}.).}
  \item{timeStrata}{
  Specified time stratification (to be chosen between \code{"year"}, \code{"quarter"}, \code{"month"}, \code{NULL},...).}
  \item{spaceStrata}{
  Specified space stratification (to be chosen between \code{"area"}, \code{"rect"}, \code{NULL},...).}
  \item{techStrata}{
  Specified technical stratification (to be chosen between \code{"gear"}, \code{"foCatNat"}, \code{"foCatEu5"}, \code{NULL},...).}
  \item{...}{Further arguments.}
}

\details{For more informations about arguments, see \emph{FishFrame/COST Exchange format specification}.}

\value{An object of class \emph{landisVol}.}

\author{Mathieu Merzereaud}

\seealso{\code{\link{landisVol}}, \code{\link{fdPlot}}, \code{\link{fdBoxplot}}, \code{\link{foPlot}}, \code{\link{foBoxplot}}
}

\examples{
data(sole)
x <- landisVolume(sole.cs,species="Solea vulgaris",fraction="LAN")
}

\keyword{methods}
