\name{fillGaps}
\alias{fillGaps}
\title{Method for completing gaps in simulated data sets}
\description{
This method completes dbeSimObj class objects with zeros in non sampled age and length dimensions 
in order to have objects of the same dimension that can be compared
}

\usage{
fillGaps(dbeSimObj, ageMin, ageMax, lenMin, lenMax)
}

\arguments{
  \item{dbeSimObj}{ A dbeOutputSim object}
  \item{ageMin}{Numeric indicating the minimum age sampled}
  \item{ageMax}{Numeric indicating the maximum age sampled}
  \item{lenMin}{Numeric indicating the minimum length sampled}
  \item{lenMax}{Numeric indicating the maximum length sampled}
}

\value{
  It returns the same object but with the age and length dimensions completed with zeros
}

\author{Dorleta Garcia \email{dgarcia@azti.es}}

%\examples{}

\keyword{methods}
