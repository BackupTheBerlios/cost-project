\name{fdPlot,landisVol-method}
\alias{fdPlot}
\alias{fdPlot,landisVol-method}
\docType{methods}
\title{Plot of mean fishing-day-catch weight by trip}
\description{
Graphical display of mean fishing-day-catch weight by trip.
It requires a \emph{landisVol} object built from \emph{landisVolume} procedure. 
Display can be done within time, technical and space stratification, if specified in \emph{landisVolume} call.
}

\usage{
fdPlot(x,groups=NULL,\dots)
}

\arguments{
  \item{x}{A \emph{landisVol} object created by \emph{landisVolume} procedure.}
  \item{groups}{
  Character specifying intra-graph stratification (to be chosen between 
  \code{"timeStrata"}, \code{"techStrata"}, \code{"spaceStrata"} and \code{NULL}).}
  \item{...}{Further graphical arguments.}
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{landisVol}}, \code{\link{landisVolume}}, \code{\link{fdBoxplot}}, \code{\link{foPlot}}, \code{\link{foBoxplot}}
}

\examples{
data(sole)

x <- landisVolume(sole.cs,species="Solea vulgaris",fraction="LAN")
fdPlot(x)

xx <- landisVolume(sole.cs,species="Solea vulgaris",fraction="LAN",
                   timeStrata="quarter",techStrata="gear")
fdPlot(xx,groups="techStrata")
}

\keyword{methods}
