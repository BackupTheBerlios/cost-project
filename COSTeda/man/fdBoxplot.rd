\name{fdBoxplot,landisVol-method}
\alias{fdBoxplot}
\alias{fdBoxplot,landisVol-method}
\docType{methods}
\title{Boxplot of raised catch weight by fishing day for each trip}
\description{
Boxplot of raised catch weight by fishing day for each trip. It requires a \emph{landisVol} object built from \emph{landisVolume} procedure.
}

\usage{
fdBoxplot(x,\dots)
}

\arguments{
  \item{x}{A \emph{landisVol} object created by \emph{landisVolume} procedure.}
  \item{...}{Further graphical arguments.}
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{landisVol}}, \code{\link{landisVolume}}, \code{\link{fdPlot}}, \code{\link{foPlot}}, \code{\link{foBoxplot}}
}

\examples{
data(sole)

x <- landisVolume(sole.cs,species="Solea vulgaris",fraction="LAN")
fdBoxplot(x)
}

\keyword{methods}
