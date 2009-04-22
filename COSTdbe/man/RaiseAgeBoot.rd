\name{RaiseAgeBoot}
\alias{RaiseAgeBoot}
\alias{RaiseAgeBoot,dbeOutput,csDataCons-method}
\docType{methods}
\title{Estimation of numbers-at-age from market sampling with bootstrap variance}
\description{
RaiseAgeBoot follows the approach of RaiseAge to calculate numbers-at-age by strata from market sampling data. A set of bootstrap estimates is created to calculate the variance.
This method requires RaiseLgthBoot to have been run to estimate the length structure and define the number of bootstrap replicates.
}

\usage{
RaiseAgeBoot(dbeOutput, csObject, type="fixed", sex=as.character(NA), bootMethod = "samples", \dots)
}
                                          
\arguments{
  \item{dbeOutput}{A \emph{dbeOutput} object containing output from RaiseLgthBoot.}
  \item{csObject}{A \emph{csDataCons} object matching 'dbeOutput' specifications.}
  \item{type}{Allocation strategy used to establish the age-length data, one of "fixed", "prop" or "ages"}
  \item{sex}{Sex}
  \item{bootMethod}{"samples" (the default) or "otoliths"}
  \item{\dots}{Further arguments}
}

\details{ If bootMethod="samples" then samples in the data are defined by PSUid and SSUid and these samples are selected with replacement to create bootstrap estimates.

 If bootMethod="otoliths" then individual age observations within each length class (i.e. row of the ALK) are resampled instead.
 This makes the assumption that individual otoliths are independent. In other words, that for age-given-length, fish in the same sample are not more similar than fish in different samples.
}                                                                               %MM

\value{An updated object of class \emph{dbeOutput}.
Slot methodDesc with bootstrap samples or bootstrap otoliths,
 nSamp\$age & nMeas\$age with number of samples and measurements,
ageStruc\$rep with bootstrap replicates for numbers-at-age, iter=0 is assigned the estimates from the original data,
 ageStruc\$estim with the mean of the bootstrap replicates,
 ageVar with the variance of the bootstrap replicates
 }



\author{David Maxwell & Mathieu Merzereaud}
\seealso{\code{\link{dbeOutput}}, \code{\link{dbeObject}}, \code{\link[COSTcore]{csDataCons}}, \code{\link[COSTcore]{clDataCons}}, \code{\link{RaiseAge}}, \code{\link{RaiseLgthBoot}} }


\examples{
data(sole)
strD <- strIni(timeStrata="quarter",techStrata="commCat")
csObject <- csDataCons(csDataVal(subset(sole.cs,sampType\%in\%c("M","V"))),strD)
clObject <- clDataCons(clDataVal(sole.cl),strD)
dbeOutput <- dbeObject(species="Solea solea",catchCat="LAN",strataDesc=strD)

# Analytical estimates
sol.dbe.an <- RaiseLgth (dbeOutput, csObject, clObject)

# Bootstrap estimates. B set to a very low number of iterations for demonstration only.
# Several thousand iterations recommended for final run, which will take some time.
sol.dbe.boot <- RaiseLgthBoot (dbeOutput, csObject, clObject, B=10)

sol.dbe.boot <- RaiseAgeBoot (dbeOutput = sol.dbe.boot, csObject = csObject, type="fixed")
}
\keyword{methods}
