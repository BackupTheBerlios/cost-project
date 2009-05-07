\name{RaiseAgebootSim}
\alias{RaiseAgeBootSim}
\docType{methods}
\title{Estimation of total numbers-at-age from market sampling with bootstrap variance for simulated data sets}
\description{
This method is the equivalent of the function \code{RaiseAgeBoot} for dbeOutputSim class objects and 
it calculates total numbers-at-age by strata from market sampling with bootstrap variance for simulated data sets
}

\usage{
RaiseAgeBootSim(dbeOutputSim,simObj,type="fixed",sex=as.character(NA),bootMethod = "samples",\dots)
}

\arguments{
  \item{dbeOutputSim}{A dbeOutputSim object.}
  \item{simObj}{A simDataCons object matching \code{dbeOutputSim} specifications.}
  \item{type}{Allocation strategy used to establish the age-length data.}
  \item{sex}{Sex}
  \item{bootMethod}{"samples" (the default) or "otoliths"}
  \item{\dots}{Further arguments}  
}

\value{An updated object of class dbeOutputSim}
Slot methodDesc with bootstrap samples or bootstrap otoliths,
 \code{nSamp\$age} & \code{nMeas\$age} with number of samples and measurements,
\code{ageStruc\$rep} with bootstrap replicates for numbers-at-age, iter=0 is assigned the estimates from the original data,
 \code{ageStruc\$estim} with the mean of the bootstrap replicates,
 \code{ageVar} with the variance of the bootstrap replicates
}

\author{Dorleta Garcia \email{dgarcia@azti.es}}

\seealso{
\code{\link{dbeOutputSim}}, \code{\link{simDataCons}}, \code{\link[COSTdbe]{RaiseAgeBoot}}
}

% \examples{ }

\keyword{methods}
