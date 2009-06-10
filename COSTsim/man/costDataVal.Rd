\name{costDataVal}
\alias{costDataVal}
\alias{costDataVal-class}
\alias{costDataVal,costData-method}
\docType{class}
\title{Class "costDataVal"}
\description{
The costDataVal-class is equivalent in structure to costData but stores validated forms of the simulated data set 
}
\section{Objects from the Class}{
The creator function "costDataVal" can be called to create objects from this class.
}

\section{Slots}{
  \tabular{lrl}{
  \bold{slot} \tab \bold{class} \tab \bold{description} \cr
  \bold{\code{desc}} \tab \code{character} \tab object description \cr
  \bold{\code{ce}} \tab \code{ceDataVal} \tab validated commercial fisheries effort data \cr
  \bold{\code{cl}} \tab \code{clDataVal} \tab validated commercial fisheries landings data \cr
  \bold{\code{cs}} \tab \code{csDataVal} \tab validated commercial fisheries sampling data \cr
  }
}

\author{Dorleta Garcia \email{dgarcia@azti.es}}

\examples{
showClass("costDataVal")
}

\keyword{classes}
