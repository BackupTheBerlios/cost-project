\name{dbePlot}
\alias{dbePlot}
\alias{dbePlot,dbeOutput-method}
\docType{methods}
\title{Graphical display of 'dbeOutput' final estimates}
\description{
Method for plotting final estimates from an input 'dbeOutput' object.
}

\usage{
dbePlot(object,Slot,type="bar",Xstratum=NULL,step=NA,dispKey=TRUE,\dots)
}

\arguments{
  \item{object}{A \emph{dbeOutput} object.}
  \item{Slot}{An 'estimates' \emph{dbeOutput} slot. If chosen value is "lenStruc", "ageStruc", "totalN" or "totalW", 'estim' element from specified slot is used.}
  \item{type}{Parameter to specify the type of the drawn plot. To be chosen between "bar" (default value), "point" and "line".}
  \item{Xstratum}{Stratum displayed on x-axis if 'Slot' is in c("nSamp", "nMes", "totalN", "totalNvar", "totalW", "totalWvar"). To be chosen between "time", "space", "technical" and \code{NULL} (default value).}
  \item{step}{Numeric. If given, empty length or age classes will be considered and displayed, according to specified value.}
  \item{dispKey}{Logical. If \code{TRUE}, a describing key is displayed}
  \item{...}{Further graphical arguments such as \emph{col, lwd, lty, pch, cex, font, rot,}\dots}
}


\author{Mathieu Merzereaud}
\seealso{\code{\link{dbeOutput}}
}

\examples{
data(sole)

#stratification object
strDef <- strIni(timeStrata="quarter",spaceStrata="area")
#consolidated object
object <- csDataCons(csDataVal(sole.cs),strDef)
#dbeOutput initial object with needed parameters
dbeOutput <- dbeObject(desc="My object",species="Solea solea",param="weight",
                       strataDesc=strDef,methodDesc="analytical")

lWeight <- bpEstim(dbeOutput,object)

dbePlot(lWeight,Slot="lenStruc",step=10,ylab="Mean weight (g)")

}
\keyword{methods}
