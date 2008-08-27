\name{NHcoast}
\alias{NHcoast}
\docType{data}
\title{ Northern Hemisphere coastline }
\description{
  1:5000000 scale costline data covereing the Northern Hemisphere. 
}
\usage{data(NHcoast)}
\format{
  A data frame with 24335 observations on the following 3 variables.
  \describe{
    \item{\code{lon}}{Longitude in decimal degrees}
    \item{\code{lat}}{Latitude in decimal degrees}
    \item{\code{string}}{Numbered sequence of coordinate strings }
  }
}
\details{
  This data set consists of sections longitude and latitude cooordinates of coastline, of variable length, separated by NA's. 
  Each section has been numbered as \$string to enable the construction of landmass polygons. 
  
}
\source{
 National Geophysical Data Center: \url{http://rimmer.ngdc.noaa.gov/coast/}
}
\seealso{\code{\link{finecoast}} as a finer resolution alternative.%modif MM 
%\code{\link{landmass.polygons}} for adding the landmass polygons.
}
\examples{
data(NHcoast)
plot(NHcoast$lon,NHcoast$lat,type="l")

%## maybe str(NHcoast) ; plot(NHcoast) ...
}
\keyword{datasets}


