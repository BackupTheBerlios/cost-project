\name{costDataCons}
\alias{costDataCons}
\alias{costDataCons-class}
\alias{costDataCons,costDataVal,strIni-method}
\alias{costDataCons,missing,missing-method}
\docType{class}
\title{Class "costDataCons"}
\description{
The costDataCons-class is equivalent in structure to costData but stores consolidated forms of the simulated data set 
}

\section{Objects from the Class}{
The creator function \code{costDataCons} can be called to create objects from this class.
}

\section{Slots}{
  \tabular{lrl}{
  \bold{slot} \tab \bold{class} \tab \bold{description} \cr
  \bold{\code{desc}} \tab \code{character} \tab object description \cr
  \bold{\code{ce}} \tab \code{ceDataCons} \tab consolidated commercial fisheries effort data \cr
  \bold{\code{cl}} \tab \code{clDataCons} \tab consolidated commercial fisheries landings data \cr
  \bold{\code{cs}} \tab \code{csDataCons} \tab consolidated commercial fisheries sampling data \cr
  }
}

\author{Dorleta Garcia \email{dgarcia@azti.es}}

\examples{
showClass("costDataCons")
}

\keyword{classes}
