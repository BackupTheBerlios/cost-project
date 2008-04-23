\name{alldepths}
\alias{alldepths}
\docType{data}
\title{Depth contours}
\description{
  50m, 100m 200m and 500m depth contours for the North Sea and Western Atlantic.   
}
\usage{data(alldepths)}
\format{
  A data frame with 2448 observations on the following 3 variables.
  \describe{
    \item{\code{lon}}{Longitude in decimal degrees}
    \item{\code{lat}}{Latitude in decimal degrees}
    \item{\code{depth}}{Depth in metres}
  }
}
\details{
Coverage is limited to (approx) -20 and 11 degrees longitude, 45 and 68 degrees latitude.
}
\source{
FRS Marine Laboratory, Aberdeen, Scotland.  
}

\examples{
data(alldepths)
data(NHcoast)
plot(NHcoast$lon,NHcoast$lat,xlim=c(-15,15),ylim=c(47,62),type="l")
lines(alldepths$lon[alldepths$depth==50],alldepths$lat[alldepths$depth==50],col=2)
lines(alldepths$lon[alldepths$depth==100],alldepths$lat[alldepths$depth==100],col=3)
lines(alldepths$lon[alldepths$depth==200],alldepths$lat[alldepths$depth==200],col=4)
lines(alldepths$lon[alldepths$depth==500],alldepths$lat[alldepths$depth==500],col=5)

%## maybe str(alldepths) ; plot(alldepths) ...
}
\keyword{datasets}

