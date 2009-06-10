\name{totVolumeSim}
\alias{totVolumeSim}
\alias{totVolumeSim,dbeOutputSim,simDataCons-method}
\docType{methods}
\title{Estimation of total volume of discards or/and landings (weight, number or number-at-length) of simulated data sets}
\description{
This function is the equivalent to \code{totVolume} for 'simDataCons' class objects. It estimates total volume of discards or/and 
landings (weight, number or number-at-length) based on various raising methods for simulated data sets.
}

\usage{
totVolumeSim(dbeOutputSim,simObj,\dots)
}

\arguments{
  \item{dbeOutputSim}{A \emph{dbeOutputSim} object. All necessary information for calculation process are taken in the first slots (species, catch category,...). See \emph{dbeObject} method for object initialization.}
  \item{simObj}{A \emph{simDataCons} object matching 'dbeOutputSim' specifications.}
  \item{...}{Further arguments such as:
  \item{type}{Specification of the raising method : \code{"trip"} (default value) for raising by trip, \code{"fo"} for raising by fishing operations, 
\code{"fd"} for raising by fishing days,\code{"landings"} for ratio-to-total landings raising method, and \code{"time"} for ratio-to-fishing duration 
raising method.}
  \item{val}{Estimated parameter. To be chosen between \code{"weight"} (default value), \code{"number"} and \code{"nAtLength"}.}
  \item{sampPar}{logical specifying if given species is considered to be automatically sampled during the sampling process (default value is \code{TRUE}).}
  \item{landSpp}{character vector describing the species considered in the 'volume of landings' variable if chosen raising method is ratio-to-landings (see 'clObject' description).}  
}
}

\value{An updated object of class dbeOutputSim.}

\references{Vigneau, J. (2006)          
\emph{Raising procedures for discards : Sampling theory (Toward agreed methodologies for calculating precision in the discard programmes)}. Working document in support of PGCCDBS (Rostock, 2006).
}

\author{Dorleta Garcia <dgarcia@azti.es>}

\seealso{
\code{\link[COSTdbe]{totVolume}}
}

%\examples{ }

\keyword{methods}
