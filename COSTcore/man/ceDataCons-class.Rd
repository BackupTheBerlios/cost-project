\name{ceDataCons-class}
\docType{class}
\alias{ceDataCons}
\alias{ceDataCons-class}
\alias{ceDataCons-methods}
\alias{ceDataCons,ceDataVal-method}
\alias{ceDataCons,missing-method}
\title{Class "ceDataCons"}
\description{This class is similar to \link{ceDataVal} but stores the information aggregated at space,time, technical level. Normaly it should be post-stratified information that can be used to apply statistical methods safely (we hope).}
\section{Objects from the Class}{
The creator function "ceDataCons" can be called to create objects from this class.
}
\section{Slots}{
	 \describe{
    \item{\code{desc}:}{Object of class \code{"character"} describing the object.}
    \item{\code{ce}:}{Object of class \code{"data.frame"} with effort information }
  }
}
\section{Methods}{
No methods defined with class "ceDataCons" in the signature.
}
\author{Ernesto Jardim <ernesto@ipimar.pt>}
\examples{
showClass("ceDataCons")
}
\keyword{classes}

