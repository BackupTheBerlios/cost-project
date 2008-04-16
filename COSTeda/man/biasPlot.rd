\name{biasPlot}
\alias{biasPlot}
\alias{biasPlot,csDataVal,clDataVal-method}
\alias{biasPlot,csDataVal,ceDataVal-method}
\alias{biasPlot,csDataCons,clDataCons-method}
\alias{biasPlot,csDataCons,ceDataCons-method}
\docType{methods}
\title{Comparative plots of Population-level/Sampling-level variables by time, technical and space strata}
\description{
This method creates an exploratory graphic to compare relative values of both a population level variable from a \emph{clDataVal/clDataCons} or a \emph{ceDataVal/ceDataCons} object, 
and a sampling level variable from a \emph{csDataVal/csDataCons} object.
Calculation can be done within time, technical and/or space stratification. 
}

\usage{
\S4method{biasPlot}{csDataVal,clDataVal}(samObj,popObj=clDataVal(),samFld="lenNum",popFld="landWt",timeStrata="quarter",spaceStrata="area",techStrata="commCat",show="all",\dots)
\S4method{biasPlot}{csDataVal,ceDataVal}(samObj,popObj=ceDataVal(),samFld="lenNum",popFld="trpNum",timeStrata="quarter",spaceStrata="area",techStrata="foCatEu5",show="all",\dots)
\S4method{biasPlot}{csDataCons,clDataCons}(samObj,popObj,samFld="lenNum",popFld="landWt",show="all",\dots)
\S4method{biasPlot}{csDataCons,ceDataCons}(samObj,popObj,samFld="lenNum",popFld="trpNum",show="all",\dots)
}

\arguments{
  \item{samObj}{A \emph{csDataVal/csDataCons} object.}
  \item{popObj}{A \emph{clDataVal/clDataCons} or a \emph{ceDataVal/ceDataCons} object.}
  \item{samFld}{A field from \emph{samObj} ("lenNum", "wt", "subSampWt" or "nbSamp" (number of samples)).}
  \item{popFld}{A field from \emph{popObj} (e.g "landWt" for 'cl' or "trpNum" for 'ce').} 
  \item{timeStrata}{
  Field specifying time stratification (e.g \code{"year"}, \code{"quarter"},\code{"month"}, \code{NULL},...). Only for 'validated' objects.}
  \item{spaceStrata}{
  Field specifying space stratification (e.g \code{"area"}, \code{"rect"}, \code{NULL},...). Only for 'validated' objects.}
  \item{techStrata}{
  Field specifying technical stratification (e.g \code{"commCat"} for \emph{clDataVal/clDataCons}, \code{"foCatNat"},
  \code{,"foCatEu5"}, \code{,"foCatEu6"}, \code{NULL},...). Only for 'validated' objects.}
  \item{show}{Character vector specifying which part to plot (to be chosen between \code{"all"}, i.e samp+pop, \code{"samp"} or \code{"pop"})}
  \item{...}{Further graphical arguments.}
}

\details{For more informations about arguments, see \emph{FishFrame/COST Exchange format specification}.}

\author{Mathieu Merzereaud}

\examples{
data(sole)
val.cs <- csDataVal(sole.cs)
val.cl <- clDataVal(sole.cl)
biasPlot(val.cs,val.cl,samFld="lenNum",popFld="landWt",timeStrata="quarter",spaceStrata="area",techStrata="commCat")
}

\keyword{methods}