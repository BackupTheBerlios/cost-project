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
FD.plot(x,groups="TechStrat",\dots)
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
data(sole3.cs)
object <- sole3.cs
#only sea sampling data is kept
object@tr <- object@tr[object@tr$sampType=="S",]
object@hh <- object@hh[object@hh$sampType=="S",]
object@sl <- object@sl[object@sl$sampType=="S",]
object@hl <- object@hl[object@hl$sampType=="S",]

x <- LD.Volume(object,fraction="LAN",species="SOL")
FD.plot(x)

xx <- LD.Volume(object,fraction="LAN",species="SOL",TimeStrat="quarter",TechStrat="gear")
FD.plot(xx)
}

\keyword{methods}