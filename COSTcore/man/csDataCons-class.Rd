\name{csDataCons-class}
\docType{class}
\alias{csDataCons}
\alias{csDataCons-class}
\alias{csDataCons-methods}
\alias{csDataCons,csDataVal-method}
\alias{csDataCons,missing-method}
\title{Class "csDataCons"}
\description{This class is similar to \link{csDataVal} but stores the information aggregated at space,time, technical level. Normaly it should be post-stratified information that can be used to apply statistical methods safely (we hope).}
\section{Objects from the Class}{
The creator function "csDataCons" can be called to create objects from this class.
}
\section{Slots}{
	 \describe{
    \item{\code{desc}:}{Object of class \code{"character"} describing the object.}
    \item{\code{tr}:}{Object of class \code{"data.frame"} with trip information.}
    \item{\code{hh}:}{Object of class \code{"data.frame"} with haul information.}
    \item{\code{sl}:}{Object of class \code{"data.frame"} with species list information.}
    \item{\code{hl}:}{Object of class \code{"data.frame"} with length information.}
    \item{\code{ca}:}{Object of class \code{"data.frame"} with growth and maturity information.}
  }
}
\section{Methods}{
No methods defined with class "csDataCons" in the signature.
}
\author{Ernesto Jardim <ernesto@ipimar.pt>}
\examples{
showClass("csDataCons")
}
\keyword{classes}

