\name{costData}
\alias{costData}
\docType{class}
\title{Class "costData"}
\description{
  The costData class stores a simulated data set
}
\section{Objects from the Class}{
The creator function \code{costData} can be called to create objects from this class.
}
\section{Slots}{
  \tabular{lrl}{
  \bold{slot} \tab \bold{class} \tab \bold{description} \cr
  \bold{\code{desc}} \tab \code{character} \tab object description \cr
  \bold{\code{ce}} \tab \code{ceData} \tab commercial fisheries effort data \cr
  \bold{\code{cl}} \tab \code{clData} \tab commercial fisheries landings data \cr
  \bold{\code{cs}} \tab \code{csData} \tab commercial fisheries sampling data \cr
  }
}
\author{Dorleta Garcia \email{dgarcia@azti.es}}

\examples{
showclass("costData")
}

\keyword{classes}
