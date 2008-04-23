\name{landmass.polygons}
\alias{landmass.polygons}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Adds landmass polygons to a map}
\description{
Adds the landmass polygons to the 1:5000000 coastline \code{NHcoast} }
\usage{
landmass.polygons(colour = "lightgrey", border = "grey", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{colour}{Colour of the landmass, default is "light grey".}
  \item{border}{Colour of the coastline, default is "grey"}
  \item{\dots}{Other arguments to polygon such as \code{density} and \code{angle}.}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\author{Alastair Pout \email{a.pout@marlab.ac.uk}}
\section{Warning }{Only suitable for use with the NGDC 1:5,000,000 scale coastline and is a bit slow. \cr Still needs some work c 28/3/2008} 
\seealso{\code{\link{NHcoast}}}
\examples{
data(NHcoast)
plot(NHcoast$lon,NHcoast$lat,xlim=c(-30,30),ylim=c(35,65),type="l")
# solid colours
landmass.polygons()
# lines
plot(NHcoast$lon,NHcoast$lat,xlim=c(-30,30),ylim=c(35,65),type="l")
landmass.polygons(density=c(5))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{aplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
