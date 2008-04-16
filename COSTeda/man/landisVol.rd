\name{landisVol-class}
\docType{class}
\alias{landisVol}
\alias{landisVol-class}
\title{Class "landis.Vol"}
\description{Sea-sampled catch informations for given species and catch category. Calculation levels are FOs, fishing days, and trips.
Any object of this class is created by \emph{landisVolume} procedure.}

\section{Slots}{
\tabular{lrl}{
\bold{slot} \tab \bold{class} \tab \bold{description} \cr
\bold{\code{fraction}} \tab \code{character} \tab Specified catch category. \cr
\bold{\code{species}} \tab \code{character} \tab Specified species. \cr

\bold{\code{Strata$timeStrata}} \tab \code{character} \tab Specified time stratification (to be chosen between \code{"year"}, \cr
                                         \tab \tab \code{"quarter"}, \code{"month"} and \code{NULL}).  \cr
\bold{\code{Strata$techStrata}} \tab \code{character} \tab Specified technical stratification (to be chosen between \code{"gear"}, \cr
                                         \tab \tab \code{"foCatNat"}, \code{"foCatEu5"}, \code{"foCatEu6"} and \code{NULL}). \cr
\bold{\code{Strata$spaceStrata}} \tab \code{character} \tab Specified technical stratification (to be chosen between \code{"area"}, \cr
                                         \tab \tab \code{"rect"} and \code{NULL}). \cr
\bold{\code{VolFO_FDTR}} \tab \code{list} \tab Catch weight by FO for each fishing day of each trip. \cr
\bold{\code{MeanFO_FDTR}} \tab \code{numeric} \tab Mean FO-catch weight for each fishing day of each trip. \cr
\bold{\code{VolFD_TR}} \tab \code{list} \tab Catch weight by fishing day for each trip and each strata, \cr
                                   \tab \tab raised by numbers of FO. \cr
\bold{\code{MeanFD_TR}} \tab \code{numeric} \tab Mean FD-raised catch weight by trip and strata.
}
}

\details{For more informations about arguments, see \emph{FishFrame/COST Exchange format specification}.}


\section{Methods}{
  \describe{
    \item{fdPlot}{\code{signature(x = "landisVol")}: graphical display of mean fishing-day-catch weight by trip (and strata if \code{groups} parameter is not \code{NULL}).}
    \item{fdBoxplot}{\code{signature(x = "landisVol")}: boxplot of catch weight by fishing day for each trip.}
    \item{foPlot}{\code{signature(x = "landisVol")}: graphical display of mean FO-catch weight for each fishing day of each trip.}
    \item{foBoxplot}{\code{signature(x = "landisVol")}: boxplot of catch weight by FO for each fishing day of each trip.}
   }
}


\author{Mathieu Merzereaud}

\seealso{\code{\link{landisVolume}}, \code{\link{fdPlot}}, \code{\link{fdBoxplot}}, \code{\link{foPlot}}, \code{\link{foBoxplot}}
}

\examples{
showClass("landisVol")
}

\keyword{classes}

