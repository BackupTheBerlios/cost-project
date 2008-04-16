\name{foPlot,landisVol-method}
\alias{foPlot}
\alias{foPlot,landisVol-method}
\docType{methods}
\title{Plot of mean sampled FO-catch weight for each fishing day of each trip}
\description{
Graphical display of mean FO-catch weight for each fishing day of each trip. It requires a \emph{landisVol} object built from \emph{landisVolume} procedure.
}

\usage{
foPlot(x,\dots)
}

\arguments{
  \item{x}{A \emph{landisVol} object created by \emph{landisVolume} procedure.}
  \item{...}{Further graphical arguments.}
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{landisVol}}, \code{\link{landisVolume}}, \code{\link{fdPlot}}, \code{\link{fdBoxplot}}, \code{\link{foBoxplot}}
}

\examples{
data(sole)

x <- landisVolume(sole.cs,species="Solea vulgaris",fraction="LAN")
foPlot(x)
}

\keyword{methods}
