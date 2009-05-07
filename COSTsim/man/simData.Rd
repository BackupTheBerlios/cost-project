\name{simData}
\alias{simData}
\docType{class}
\title{Class "simData"}
\description{
  The simData class stores various costData class objects simulated according to the Bayesian model in library \code{COSTmbe}
}

\section{Objects from the Class}{
The creator function \code{simData} can be called to create objects from this class.
}

\section{Slots}{
  \tabular{lrl}{
  \bold{slot} \tab \bold{class} \tab \bold{description} \cr
  \bold{\code{desc}} \tab \code{character} \tab object description \cr
  \bold{\code{species}} \tab \code{character} \tab species description. Recall of SL\$spp and SL\$sex \cr
  \bold{\code{samples}} \tab \code{list} \tab each element in the list is a costData object representing a simulated data set \cr
  \bold{\code{initial.fit}} \tab \code{list} \tab  \cr
  \bold{\code{setup.args}} \tab \code{list} \tab set up parameters \cr
  \bold{\code{burnin}} \tab \code{numeric} \tab  \cr
  \bold{\code{nmcmc}} \tab \code{numeric} \tab  \cr
  \bold{\code{l.int}} \tab \code{numeric} \tab  \cr
  \bold{\code{Int}} \tab \code{list} \tab  \cr
  \bold{\code{Slp}} \tab \code{list} \tab  \cr
  \bold{\code{landings}} \tab \code{numeric} \tab  \cr
  \bold{\code{nHaul}} \tab \code{integer} \tab  \cr
  \bold{\code{nseas}} \tab \code{numeric} \tab  \cr
  }
}

\author{Dorleta García \email{dgarcia@azti.es}}

\examples{
showclass("simData")
}

\keyword{classes}
