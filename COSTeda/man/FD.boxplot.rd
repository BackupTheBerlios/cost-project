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
data(sole)
object <- sole.cs

x <- LD.Volume(object,fraction="LAN",species="Solea vulgaris")
FD.boxplot(x)
}

\keyword{methods}
