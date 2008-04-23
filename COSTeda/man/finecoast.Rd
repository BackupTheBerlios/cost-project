\name{finecoast}
\alias{finecoast}
\docType{data}
\title{ Western European coastline}
\description{
  1:250,000 scale coastline data for Western Europe.  
  \cr Extends to longitude 15 degrees east and about -14 west, latitude 45 to 63, hence mainly covers UK and North Sea
}
\usage{data(finecoast)}
\format{
  A data frame with 429912 observations on the following 2 variables.
  \describe{
    \item{\code{lon}}{Longitude in decimal degrees}
    \item{\code{lat}}{Latitude in decimal degrees}
  }
}
%\details{
%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
 National Geophysical Data Center \url{http://rimmer.ngdc.noaa.gov/coast/}
}
\seealso{\code{\link{NHcoast}}}
%\references{
 % ~~ possibly secondary sources and usages ~~
%}
\examples{
data(finecoast)
plot(finecoast$lon,finecoast$lat,type="l")

%## maybe str(finecoast) ; plot(finecoast) ...
}
\keyword{datasets}
