\name{dbePlot}
\alias{dbePlot}
\alias{dbePlot,dbeOutput-method}
\docType{methods}
\title{Graphical display of 'dbeOutput' final estimates}
\description{
Method for plotting final estimates from an input 'dbeOutput' object.
}

\usage{
dbePlot(object,Slot,type="bar",Xstratum=NULL,dispKey=TRUE,\dots)
}

\arguments{
  \item{object}{A \emph{dbeOutput} object.}
  \item{Slot}{An 'estimates' \emph{dbeOutput} slot. If chosen value is "lenStruc", "ageStruc", "totalN" or "totalW", 'estim' element from specified slot is used.}
  \item{type}{Parameter to specify the type of the drawn plot. To be chosen between "bar" (default value), "point" and "line".}
  \item{Xstratum}{Stratum displayed on x-axis if 'Slot' is in c("nSamp","nMes","totalN","totalNvar","totalW","totalWvar"). To be chosen between "time", "space", "technical" and \code{NULL} (default value).}
  \item{dispKey}{Logical. If \code{TRUE}, a describing key is displayed}
  \item{...}{Further graphical arguments such as \emph{col, lwd, lty, pch, cex, font, rot,}\dots}
}


\author{Mathieu Merzereaud}
\seealso{\code{\link{dbeOutput}}
}

\examples{
data(sole)

#consolidated datasets are built 
strDef <- strIni(timeStrata="quarter",techStrata="foCatEu5")
csObject <- csDataCons(csDataVal(sole.cs),strDef)
clObject <- clDataCons(clDataVal(sole.cl),strDef)
ceObject <- ceDataCons(ceDataVal(sole.ce),strDef)

#dbeOutput initial object
obj <- dbeObject(desc="My object",species="Solea solea",catchCat="DIS",strataDesc=strDef,methodDesc="analytical")

#raising by trip
newObj <- totVolume(obj,csObject,ceObject)

#plotting method
dbePlot(newObj,"totalW",type="bar",Xstratum="technical")

}
\keyword{methods}
