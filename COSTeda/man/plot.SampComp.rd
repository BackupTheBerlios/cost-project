\name{SampComp.plot,csData,clData-method}
\alias{SampComp.plot}
\alias{SampComp.plot,csData,clData-method}
\docType{methods}
\title{Plots of volume of landings and number of fish measured by time, technical and space strata}
\description{
This method creates an exploratory graphic to compare relative volume of landings and relative number of fish measured.
Calculation can be done for time, technical and/or space stratification. It requires a \emph{csData} object and a \emph{clData} object built from \pkg{COSTcore} package.
}

\usage{
SampComp.plot(object1,object2,TimeStrat="quarter",TechStrat="commCat",
              SpaceStrat="area",show="all",separate=FALSE,\dots)
}

\arguments{
  \item{object1}{A \emph{csData} object with \emph{tr}, \emph{hh}, \emph{sl} and \emph{hl} informations.}
  \item{object2}{A \emph{clData} object. 'cl' table data defines stratas occurrences.}
  \item{TimeStrat}{
  Field specifying time stratification (to be chosen between \code{"year"}, \code{"quarter"},
  \code{"month"} and \code{NULL}). See Details.}
  \item{TechStrat}{
  Field specifying technical stratification (to be chosen between \code{"commCat"}, \code{"foCatNat"},
  \code{,"foCatEu5"}, \code{,"foCatEu6"} and \code{NULL}). See Details.}
  \item{SpaceStrat}{
  Field specifying space stratification (to be chosen between \code{"area"}, \code{"rect"}
  and \code{NULL}). See Details.}
  \item{show}{Character vector specifying which stratas to plot (to be chosen between \code{"Time"}, \code{"Technical"},
  \code{"Space"} and/or \code{"all"})}
  \item{separate}{
  Logical value to specify the treatment of missing values in each strata. \code{FALSE} means that
  calculation will be done within all stratas. \code{TRUE} means that it will be done within each strata separately.}
  \item{...}{Further graphical arguments.}
}

\details{For more informations about arguments, see \emph{FishFrame/COST Exchange format specification}.}

\author{Mathieu Merzereaud}

\examples{
data(sole3.cs)
data(sole3.cl)
SampComp.plot(sole3.cs,sole3.cl)
}

\keyword{methods}