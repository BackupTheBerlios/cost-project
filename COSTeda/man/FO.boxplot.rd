\name{FO.boxplot,LD.Vol-method}
\alias{FO.boxplot}
\alias{FO.boxplot,LD.Vol-method}
\docType{methods}
\title{Boxplot of sampled catch weight by FO for each fishing day of each trip}
\description{
Boxplot of catch weight by FO for each fishing day of each trip. It requires a \emph{LD.Vol} object built from \emph{LD.Volume} procedure.
}

\usage{
FO.boxplot(x,\dots)
}

\arguments{
  \item{x}{A \emph{LD.Vol} object created by \emph{LD.Volume} procedure.}
  \item{...}{Further graphical arguments.}
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{LD.Vol}}, \code{\link{LD.Volume}}, \code{\link{FD.plot}}, \code{\link{FD.boxplot}}, \code{\link{FO.plot}}
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
FO.boxplot(x)
}

\keyword{methods}
