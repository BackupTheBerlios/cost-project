\name{RaiseAgeSim}
\alias{RaiseAgeSim}
\alias{RaiseAgeSim-methods}
\alias{RaiseAgeSim,dbeOutputSim,simDataCons-method}
\docType{methods}
\title{Estimation of total numbers-at-age from market sampling for simulated data sets}
\description{
This method is the equivalent of the function \code{RaiseAge} for dbeOutputSim class objects 
and it calculates total numbers-at-age by strata from market sampling for simulated data sets
}

\usage{
RaiseAgeSim(dbeOutputSim,simObj,type="p",sex=as.character(NA),\dots)
}

\arguments{
  \item{dbeOutputSim}{A dbeOutputSim object.}
  \item{simObj}{A simDataCons object matching \code{dbeOutputSim} specifications.}
  \item{type}{Allocation strategy used to establish the age-length data.}
  \item{sex}{Sex}
  \item{\dots}{Further arguments}  
}

\value{An updated object of class dbeOutputSim.
Slots \code{nSamp\$age} & \code{nMeas\$age} with number of samples and measurements,
\code{ageStruc\$estim} with numbers-at-age estimates,
\code{ageVar} with the variance of numbers-at-age .
 }

\author{Dorleta Garcia \email{dgarcia@azti.es}}

\seealso{
\code{\link{dbeOutputSim}}, \code{\link{simDataCons}}, \code{\link[COSTdbe]{RaiseAge}}
}

%\examples{}

\keyword{methods}
