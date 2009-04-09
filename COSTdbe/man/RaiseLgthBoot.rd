\name{RaiseLgthBoot}
\alias{RaiseLgthBoot}
\alias{RaiseLgth,dbeOutput,csDataCons,clDataCons-method}
\docType{methods}
\title{Estimation of numbers-at-length from market sampling with bootstrap variance}
\description{
RaiseLgthBoot follows the approach of RaiseLgth to calculate numbers-at-length by strata from market sampling data. The data are resampled with replacement to create B bootstrap estimates that are used to calculate the variance.
}

\usage{
RaiseLgthBoot(dbeOutput,csObject,clObject,spp,taxon,sex=as.character(NA),B ,\dots)
}

\arguments{
  \item{dbeOutput}{A \emph{dbeOutput} object.}
  \item{csObject}{A \emph{csDataCons} object matching 'dbeOutput' specifications.}
  \item{clObject}{A \emph{clDataCons} object matching 'dbeOutput' specifications.}
  \item{spp}{Species, if missing this is set to dbeOutput@species}
  \item{taxon}{Taxon, if missing this is set to dbeOutput@species}
  \item{sex}{Sex}
  \item{B}{Number of bootstrap interations}
  \item{\dots}{Further arguments}  
}

\references{}

\value{An updated object of class \emph{dbeOutput}.
Slots nSamp$len & nMeas$len with number of samples and measurements,
 methodDesc with "bootstrap",
 totalW$estim with total weight,
 lenStruc$rep & totalN$rep with bootstrap replicates for numbers-at-length and total numbers, iter=0 is assigned the estimates from the original data,
 lenStruc$estim & totalN$estim with the mean of the bootstrap replicates,
 lenVar & totalNvar with the variance of the bootstrap replicates}


\author{David Maxwell & Mathieu Merzereaud}
\seealso{\code{\link{dbeOutput}}, \code{\link{dbeObject}}, \code{\link[COSTcore]{csDataCons}}, \code{\link[COSTcore]{clDataCons}}, \code{\link{RaiseLgth}} }
}

\examples{
data(sole)
strD <- strIni(timeStrata="quarter",techStrata="commCat")
csObject <- csDataCons(csDataVal(subset(sole.cs,sampType%in%c("M","V"))),strD)
clObject <- clDataCons(clDataVal(sole.cl),strD)
dbeOutput <- dbeObject(species="Solea solea",catchCat="LAN",strataDesc=strD)

# Analytical estimates
sol.dbe.an <- RaiseLgth (dbeOutput, csObject, clObject)

# Bootstrap estimates. B set to a very low number of iterations for demonstration only.
# Several thousand iterations recommended for final run, which will take some time.
sol.dbe.boot <- RaiseLgthBoot (dbeOutput, csObject, clObject, B=10)

}
\keyword{methods}
