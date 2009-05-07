\name{RaiseLgthSim}
\alias{RaiseLgthSim}
\docType{methods}
\title{Estimation of total numbers-at-length from market sampling for simulated data sets}
\description{
This method is the equivalent of the function \code{RaiseLgth} for dbeOutputSim class objects 
and it calculates total numbers-at-length by strata from market sampling for simulated data sets
}

\usage{
RaiseLgthSim(dbeOutputSim,simObj,spp,taxon,sex=as.character(NA),\dots)
}

\arguments{
  \item{dbeOutputSim}{A dbeOutputSim object.}
  \item{simObj}{A simDataCons object matching 'dbeOutputSim' specifications.}
  \item{spp}{Species, if missing this is set to dbeOutput@species}
  \item{taxon}{Taxon, if missing this is set to dbeOutput@species}
  \item{sex}{Sex}
  \item{\dots}{Further arguments}  
}

\value{An updated object of class dbeOutputSim.
Slots \code{nSamp\$len} & \code{nMeas\$len} with number of samples and measurements,
\code{totalW\$estim} with total weight,
\code{lenStruc\$estim} with numbers-at-length estimates,
\code{lenVar} with the variance of numbers-at-length.
 }

\author{Dorleta Garcia <dgarcia@azti.es>}

\seealso{
\code{\link{dbeOutputSim}}, \code{\link{simDataCons}}, \code{\link[COSTdbe]{RaiseLgth}}
}

%\examples{}

\keyword{methods}
