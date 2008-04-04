\name{SampComp.plot}
\alias{SampComp.plot}
\alias{SampComp.plot,csData,clData-method}
\alias{SampComp.plot,csData,ceData-method}
\docType{methods}
\title{Comparative plots of Pop/Sampled variables by time, technical and space strata}
\description{
This method creates an exploratory graphic to compare relative values of both a population level variable from a \emph{clData} or a \emph{ceData} object, 
and a sampling level variable from a \emph{csData} object.
Calculation can be done within time, technical and/or space stratification (if they're available in both population and sampling tables). 
It requires a \emph{csData} object and a \emph{clData} object built from \pkg{COSTcore} package.
}

\usage{
SampComp.plot(object1,object2,var1="lenNum",var2,TimeStrat="quarter",TechStrat,
              SpaceStrat="area",show="all",\dots)
}

\arguments{
  \item{object1}{A \emph{csData} object with \emph{tr}, \emph{hh}, \emph{sl} and \emph{hl} informations.}
  \item{object2}{A \emph{clData} or a \emph{ceData} object.}
  \item{var1}{A numerical field from \emph{object1}.}
  \item{var2}{A numerical field from \emph{object2}.} 
  \item{TimeStrat}{
  Field specifying time stratification (e.g \code{"year"}, \code{"quarter"},
  \code{"month"}, \code{NULL},...). See Details.}
  \item{TechStrat}{
  Field specifying technical stratification (e.g \code{"commCat"} for \emph{clData}, \code{"foCatNat"},
  \code{,"foCatEu5"}, \code{,"foCatEu6"}, \code{NULL},...). See Details.}
  \item{SpaceStrat}{
  Field specifying space stratification (e.g \code{"area"}, \code{"rect"}, \code{NULL},...). See Details.}
  \item{show}{Character vector specifying which part to plot (to be chosen between \code{"all"}, e.g samp+pop, \code{"samp"} or \code{"pop"})}
  \item{...}{Further graphical arguments.}
}

\details{For more informations about arguments, see \emph{FishFrame/COST Exchange format specification}.}

\author{Mathieu Merzereaud}

\examples{
data(sole)
SampComp.plot(sole.cs,sole.cl,var1="lenNum",var2="landWt",TimeStrat="quarter",TechStrat="commCat",SpaceStrat="area",show="all")
}

\keyword{methods}