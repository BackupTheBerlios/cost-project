\name{bpEstim}
\alias{bpEstim}
\alias{bpEstim,dbeOutput,csDataCons-method}
\docType{methods}
\title{Analytical estimates of biological parameters}
\description{
This method implements analytical estimates of empirical weight-at-length/age, maturity-at-length/age, sex-ratio-at-length/age and variances. The needed parameters are from input 'dbeOutput' slots. 
}                                                                                                                                           

\usage{
bpEstim(dbeOutput,object,adjust=TRUE,immature.scale=1,\dots)
}

\arguments{
  \item{dbeOutput}{A \emph{dbeOutput} object.}
  \item{object}{A \emph{csDataCons} object.}
  \item{adjust}{Logical. If FALSE, length distribution in \emph{object}'s CA table is supposed to be representative of the catch (all calculations are made within CA). 
If TRUE (default value), previous assumption is rejected, and estimates-at-age are calculated by injecting \emph{object}'s HL information.}
  \item{immature.scale}{Numeric or character. Specifies the value(s) in \emph{matStage} field (from ca table in \emph{object}) for which the individuals are defined as immature.}
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

}
\keyword{methods}
