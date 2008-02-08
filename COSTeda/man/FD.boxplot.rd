\name{FD.boxplot,LD.Vol-method}
\alias{FD.boxplot}
\alias{FD.boxplot,LD.Vol-method}
\docType{methods}
\title{Boxplot of raised catch weight by fishing day for each trip}
\description{
Boxplot of raised catch weight by fishing day for each trip. It requires a \emph{LD.Vol} object built from \emph{LD.Volume} procedure.
}

\usage{
FD.boxplot(x,\dots)
}

\arguments{
  \item{x}{A \emph{LD.Vol} object created by \emph{LD.Volume} procedure.}
  \item{...}{Further graphical arguments.}
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{LD.Vol}}, \code{\link{LD.Volume}}, \code{\link{FD.plot}}, \code{\link{FO.plot}}, \code{\link{FO.boxplot}}
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
FD.boxplot(x)
}

\keyword{methods}