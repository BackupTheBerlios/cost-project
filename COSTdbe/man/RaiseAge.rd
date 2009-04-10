\name{RaiseAge}
\alias{RaiseAge}
\alias{RaiseAge,dbeOutput,csDataCons-method}
\docType{methods}
\title{Estimation of total numbers-at-age from market sampling}
\description{
This method calculates total numbers-at-age by strata from market sampling data. 
As calculation requires total numbers-at-length estimates, \emph{RaiseLgth} method must be run beforehand to update input \emph{dbeOutput} object with needed information.  
}

\usage{
RaiseAge(dbeOutput,csObject,type="fixed",sex=as.character(NA),\dots)
}

\arguments{
  \item{dbeOutput}{A \emph{dbeOutput} object.}
  \item{csObject}{A \emph{csDataCons} object matching 'dbeOutput' specifications.}
  \item{type}{Allocation strategy used to establish the age-length data.}
  \item{sex}{Sex}
  \item{\dots}{Further arguments}  
}

%\references{}

\value{An updated object of class \emph{dbeOutput}.
Slots nSamp\$age & nMeas\$age with number of samples and measurements,
ageStruc\$estim with numbers-at-age estimates,
ageVar with the variance of numbers-at-age .
 }


\author{Mathieu Merzereaud}
\seealso{\code{\link{dbeOutput}}, \code{\link{dbeObject}}, \code{\link{RaiseLgth}}, \code{\link[COSTcore]{csDataCons}}, \code{\link[COSTcore]{clDataCons}}
}

\examples{
data(sole)
#stratification
strD <- strIni(timeStrata="quarter",techStrata="commCat")
#only market sampling data and biological parameters are kept
csObject <- csDataCons(csDataVal(subset(sole.cs,sampType\%in\%c("M","V"))),strD)
clObject <- clDataCons(clDataVal(sole.cl),strD)
#initializing the output object
dbeOutput <- dbeObject(species="Solea solea",catchCat="LAN",strataDesc=strD)

# total numbers at length
dbeOutput <- RaiseLgth (dbeOutput, csObject, clObject)

# total numbers at age
dbeOutput <- RaiseAge (dbeOutput, csObject)


}
\keyword{methods}
