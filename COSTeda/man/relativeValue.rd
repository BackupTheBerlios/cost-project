\name{relativeValue}
\alias{relativeValue}
\alias{relativeValue,csDataVal,missing-method}
\alias{relativeValue,csDataVal,strIni-method}
\alias{relativeValue,clDataVal,missing-method}
\alias{relativeValue,clDataVal,strIni-method}
\alias{relativeValue,ceDataVal,missing-method}
\alias{relativeValue,ceDataVal,strIni-method}
\alias{relativeValue,csDataCons,missing-method}
\alias{relativeValue,clDataCons,missing-method}
\alias{relativeValue,ceDataCons,missing-method}
\docType{methods}
\title{Calculation of relative values of a numerical field within time, space and/or technical stratification}
\description{
This method calculates relative values of a population level variable (if input object class is \emph{clDataVal/clDataCons/ceDataVal/ceDataCons}),
or a sampling level variable (if input object class is \emph{csDataVal/csDataCons}).
Calculation can be done within time, space and/or technical stratification. 
Output is an \emph{edaResult} object with \emph{desc=}"csRelativeValue" or "clceRelativeValue".
An exploratory graphic to compare two objects can be made by applying \emph{plot} function. (see plot.edaResult)
}

\usage{
\S4method{relativeValue}{csDataVal,missing}(data,field="lenNum",\dots)
\S4method{relativeValue}{csDataVal,strIni}(data,strDef,field="lenNum",\dots)
\S4method{relativeValue}{clDataVal,missing}(data,field="landWt",\dots)
\S4method{relativeValue}{clDataVal,strIni}(data,strDef,field="landWt",\dots)
\S4method{relativeValue}{csDataVal,missing}(data,field="trpNum",\dots)
\S4method{relativeValue}{csDataVal,strIni}(data,strDef,field="trpNum",\dots)
\S4method{relativeValue}{csDataCons,missing}(data,field="lenNum",\dots)
\S4method{relativeValue}{clDataCons,missing}(data,field="landWt",\dots)
\S4method{relativeValue}{ceDataCons,missing}(data,field="trpNum",\dots)
}

\arguments{
  \item{data}{An \emph{csDataVal/clDataVal/ceDataVal/csDataCons/clDatacons/ceDataCons} object.}
  \item{strDef}{Optionnal. A \emph{strIni} object. Specified stratification must match with 'data' parameter.}
  \item{field}{A numerical field from \emph{data} (e.g "lenNum", "wt", "subSampWt" or "nbSamp" (number of samples) for 'cs', 
  "landWt" for 'cl', or "trpNum" for 'ce').}
  \item{...}{Further arguments.}
}


\author{Mathieu Merzereaud}

\examples{
data(sole)
sole.cs.val <- csDataVal(sole.cs)
sole.cl.val <- clDataVal(sole.cl) 
strD <- strIni(timeStrata="month",spaceStrata="area",techStrata="commCat")

CS <- relativeValue(sole.cs.val,strD,"nbSamp")
CL <- relativeValue(sole.cl.val,strD)
}

\keyword{methods}