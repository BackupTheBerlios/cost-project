\name{foBoxplot,landisVol-method}
\alias{foBoxplot}
\alias{foBoxplot,landisVol-method}
\docType{methods}
\title{Boxplot of sampled catch weight by FO for each fishing day of each trip}
\description{
Boxplot of catch weight by FO for each fishing day of each trip. It requires a \emph{landisVol} object built from \emph{landisVolume} procedure.
}

\usage{
foBoxplot(x,\dots)
}

\arguments{
  \item{x}{A \emph{landisVol} object created by \emph{landisVolume} procedure.}
  \item{...}{Further graphical arguments.}
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{landisVol}}, \code{\link{landisVolume}}, \code{\link{fdPlot}}, \code{\link{fdBoxplot}}, \code{\link{foPlot}}
}

\examples{
data(sole)

x <- landisVolume(sole.cs,species="Solea vulgaris",fraction="LAN")
foBoxplot(x)
}

\keyword{methods}
