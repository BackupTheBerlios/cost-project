\name{bpEstim}
\alias{bpEstim}
\alias{bpEstim,dbeOutput,csDataCons-method}
\docType{methods}
\title{Analytical estimates of biological parameters}
\description{
This method implements analytical estimates of empirical weight-at-length, maturity-at-length, sex-ratio-at-length and variances. The needed parameters are from input 'dbeOutput' slots. 
}                                                                                                                                           

\usage{
bpEstim(dbeOutput,object,\dots)
}

\arguments{
  \item{dbeOutput}{A \emph{dbeOutput} object.}
  \item{object}{A \emph{csDataCons} object.}
  \item{...}{Further arguments.}
}


\author{Mathieu Merzereaud}
\seealso{\code{\link{dbeOutput},\link[COSTcore]{csDataCons}}
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
#graphic
dbePlot(lWeight,Slot="lenStruc",step=10,ylab="Mean weight (g)")

}
\keyword{methods}
