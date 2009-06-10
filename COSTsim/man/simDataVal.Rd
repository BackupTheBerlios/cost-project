\name{simDataVal}
\alias{simDataVal}
\alias{simDataVal-class}
\alias{simDataVal,simData-method}
\docType{class}
\title{Class "simDataVal"}
\description{
  The simDataVal-class is equivalent in structure to simData but stores validated forms of the simulated data sets 
}

\section{Objects from the Class}{
The creator function \code{simDataVal} can be called to create objects from this class.
}

\section{Slots}{
  \tabular{lrl}{
  \bold{slot} \tab \bold{class} \tab \bold{description} \cr
  \bold{\code{desc}} \tab \code{character} \tab object description \cr
  \bold{\code{species}} \tab \code{character} \tab species description. Recall of SL\$spp and SL\$sex \cr
  \bold{\code{samples}} \tab \code{list} \tab each element in the list is a costDataVal object representing a simulated data set \cr
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

\author{Dorleta Garcia \email{dgarcia@azti.es}}

\examples{
showClass("simDataVal")
}

\keyword{classes}
