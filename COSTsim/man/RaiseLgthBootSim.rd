\name{RaiseLgthBootSim}
\alias{RaiseLgthBootSim}
\docType{methods}
\title{Estimation of total numbers-at-length from market sampling with bootstrap variance for simulated data sets}
\description{
This method is the equivalent of the function \code{RaiseLgthBoot} for dbeOutputSim class objects and 
it calculates total numbers-at-length by strata from market sampling with bootstrap variance for simulated data sets
}

\usage{
RaiseLgthBootSim(dbeOutputSim,simObj,spp,taxon,sex=as.character(NA),B,\dots)
}

\arguments{
  \item{dbeOutputSim}{A dbeOutputSim object.}
  \item{simObj}{A simDataCons object matching \code{dbeOutputSim} specifications.}
  \item{spp}{Species, if missing this is set to dbeOutput@species}
  \item{taxon}{Taxon, if missing this is set to dbeOutput@species}
  \item{sex}{Sex}
  \item{B}{Number of bootstrap interations}
  \item{\dots}{Further arguments}  
}


\value{An updated object of class dbeOutputSim.
Slots \code{nSamp\$len} & \code{nMeas\$len} with number of samples and measurements,
 \code{methodDesc} with "bootstrap",
 \code{totalW\$estim} with total weight,
 \code{lenStruc\$rep} & \code{totalN\$rep} with bootstrap replicates for numbers-at-length and total numbers, iter=0 is assigned the estimates from the original data,
 \code{lenStruc\$estim} & \code{totalN\$estim} with the mean of the bootstrap replicates,
 \code{lenVar} & \code{totalNvar} with the variance of the bootstrap replicates
 }

\author{Dorleta Garcia <dgarcia@azti.es>}

\seealso{
\code{\link{dbeOutputSim}}, \code{\link{simDataCons}}, \code{\link[COSTdbe]{RaiseLgthBoot}}
}

%\examples{ }

\keyword{methods}
