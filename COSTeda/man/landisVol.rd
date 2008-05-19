\name{landisVol}
\alias{landisVol}
\alias{landisVol,csData,missing-method}
\alias{landisVol,csData,strIni-method}
\docType{methods}
\title{Calculation upon sea-sampled catch at FO level, fishing day level, and trip level.}
\description{
This method creates an object of class \emph{edaResult} with \emph{desc="landisVol"} containing volume informations about sampled and raised catch,
for a given species in a specified catch category. It requires a \emph{csData} object built from \pkg{COSTcore} package.
Only sea sampling data is computed.
}


\usage{
\S4method{landisVol}{csData,missing}(object,species,fraction="LAN",\dots)
\S4method{landisVol}{csData,strIni}(object,strDef,species,fraction="LAN",\dots)
}

\arguments{
  \item{object}{A \emph{csData} object with sea-sampling information (\emph{tr}, \emph{hh} and \emph{sl} required).}
  \item{strDef}{A \emph{strIni} object specifying time (e.g \code{"year"}, \code{"quarter"}, \code{"month"},...),
  space (e.g \code{"area"}, \code{"rect"},...) and/or technical stratification  (e.g \code{"gear"}, \code{"foCatNat"}, \code{"foCatEu5"},...).}
  \item{species}{Field specifying species (e.g \code{"Solea vulgaris"}).}
  \item{fraction}{Field specifying catch category (to be chosen between \code{"LAN"} and \code{"DIS"}.).}
  \item{...}{Further arguments.}
}


\value{An object of class \emph{edaResult} with \emph{desc="landisVol"}.}

\author{Mathieu Merzereaud}

\seealso{\code{\link{edaResult}}, \code{\link{plot.edaResult}}, \code{\link{boxplot.edaResult}}
}

\examples{
data(sole)  
obj <- landisVol(sole.cs,species="Solea vulgaris")
}

\keyword{methods}
