\name{SampComp.plot,csData,clData-method}
\alias{SampComp.plot}
\alias{SampComp.plot,csData,clData-method}
\docType{methods}
\title{Comparative plots of Pop/Sampled variables by time, technical and space strata}
\description{
This method creates an exploratory graphic to compare relative values of both a population level variable from a \emph{clData} object, 
and a sampling level variable from a \emph{csData} object.
Calculation can be done within time, technical and/or space stratification. It requires a \emph{csData} object and a \emph{clData} object built from \pkg{COSTcore} package.
}

\usage{
SampComp.plot(object1,object2,var1="lenNum",var2="landWt",TimeStrat="quarter",TechStrat="commCat",
              SpaceStrat="area",show="all",\dots)
}

\arguments{
  \item{object1}{A \emph{csData} object with \emph{tr}, \emph{hh}, \emph{sl} and \emph{hl} informations.}
  \item{object2}{A \emph{clData} object. 'cl' table data defines stratas occurrences.}
  \item{var1}{A numerical field from \emph{csData} object.}
  \item{var2}{A \emph{clData} object.} 
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
  \item{...}{Further graphical arguments.}
}

\details{For more informations about arguments, see \emph{FishFrame/COST Exchange format specification}.}

\author{Mathieu Merzereaud}

\examples{
data(sole)
SampComp.plot(sole.cs,sole.cl,var1="lenNum",var2="landWt",TimeStrat="quarter",TechStrat="commCat",SpaceStrat="area",show="all")
}

\keyword{methods}