\name{FD.plot,LD.Vol-method}
\alias{FD.plot}
\alias{FD.plot,LD.Vol-method}
\docType{methods}
\title{Plot of mean fishing-day-catch weight by trip}
\description{
Graphical display of mean fishing-day-catch weight by trip.
It requires a \emph{LD.Vol} object built from \emph{LD.Volume} procedure. 
Display can be done within time, technical and space stratification, if specified in \emph{LD.Volume} call.
}

\usage{
FD.plot(x,groups=NULL,\dots)
}

\arguments{
  \item{x}{A \emph{LD.Vol} object created by \emph{LD.Volume} procedure.}
  \item{groups}{
  Character specifying intra-graph stratification (to be chosen between 
  \code{"TimeStrat"}, \code{"TechStrat"}, \code{"SpaceStrat"} and \code{NULL}).}
  \item{...}{Further graphical arguments.}
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{LD.Vol}}, \code{\link{LD.Volume}}, \code{\link{FD.boxplot}}, \code{\link{FO.plot}}, \code{\link{FO.boxplot}}
}

\examples{
data(sole)
object <- sole.cs

x <- LD.Volume(object,fraction="LAN",species="Solea vulgaris")
FD.plot(x)

xx <- LD.Volume(object,fraction="LAN",species="Solea vulgaris",TimeStrat="quarter",TechStrat="gear")
FD.plot(xx,groups="TechStrat")
}

\keyword{methods}
