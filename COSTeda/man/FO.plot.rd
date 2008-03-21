\name{FO.plot,LD.Vol-method}
\alias{FO.plot}
\alias{FO.plot,LD.Vol-method}
\docType{methods}
\title{Plot of mean sampled FO-catch weight for each fishing day of each trip}
\description{
Graphical display of mean FO-catch weight for each fishing day of each trip. It requires a \emph{LD.Vol} object built from \emph{LD.Volume} procedure.
}

\usage{
FO.plot(x,\dots)
}

\arguments{
  \item{x}{A \emph{LD.Vol} object created by \emph{LD.Volume} procedure.}
  \item{...}{Further graphical arguments.}
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{LD.Vol}}, \code{\link{LD.Volume}}, \code{\link{FD.plot}}, \code{\link{FD.boxplot}}, \code{\link{FO.boxplot}}
}

\examples{
data(sole)
object <- sole.cs
#only sea sampling data is kept
object@tr <- object@tr[object@tr$sampType=="S",]
object@hh <- object@hh[object@hh$sampType=="S",]
object@sl <- object@sl[object@sl$sampType=="S",]
object@hl <- object@hl[object@hl$sampType=="S",]

x <- LD.Volume(object,fraction="LAN",species="Solea vulgaris")
FO.plot(x)
}

\keyword{methods}
