\name{ICESAreaRects}
\alias{ICESAreaRects}
\docType{data}
\title{ICES areas and statistical rectangles}
\description{
Data frame containing ICES area codes and ICES statistical rectangle codes.
}
\usage{data(ICESAreaRects)}
\format{
  A data frame with 4403 observations on the following 8 variables:
  \describe{
  \item{\code{TypeName}}{area type (all ICES)}
  \item{\code{AreaCode}}{area code (to lowest level hence same as subdivision)}
  \item{\code{AreaName}}{area name}
  \item{\code{StatRect}}{statistical rectangle (including some subdivisions e.g. "51E2a")}
  %added MM
  \item{\code{lon}}{longitude}
  \item{\code{lat}}{latitude}
  %added MM
  \item{\code{statsq}}{statistical rectangle (without subdivisions)}
  \item{\code{division}}{ICES division}
  \item{\code{subarea}}{ICES subarea}  
  \item{\code{subdivision}}{ICES subdivision}  
}   }
%\details{
 % ~~ If necessary, more details than the __description__ above ~~
%}
\source{
  FRS Marine Laboratory, Aberdeen, Scotland. 
  }
\examples{
data(ICESAreaRects)
head(ICESAreaRects)

}
\keyword{datasets}
