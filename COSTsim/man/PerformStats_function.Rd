\name{PerformStats-methods}
\alias{PerformStats-methods}
\alias{PerformStats,dbeOutputSim,trueData-method}
\alias{PerformStats,missing,missing-method}
\docType{methods}
\title{Method for computing performance statistics}
\description{
  The function computes the performance statistics to evaluate accuracy, bias and precision of different estimators. The outcome 
  object is of the performStats class. 
}
\usage{
performStats(estSimObj, trueDataObj, desc, nSamples)
}
\arguments{
  \item{estSimObj}{A dbeOutputSim object}
  \item{trueDataObj}{A trueData object}
  \item{desc}{Object description}
  \item{nSamples}{Not needed anymore?}
}

\value{
An object of the class performStats
}

\author{Dorleta Garcia <dgarcia@azti.es>}

\seealso{
\code{\link{PerformStats}}
}

%\examples{}

\keyword{methods}
