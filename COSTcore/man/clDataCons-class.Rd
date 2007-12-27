\name{clDataCons-class}
\docType{class}
\alias{clDataCons}
\alias{clDataCons-class}
\alias{clDataCons-methods}
\alias{clDataCons,clDataVal-method}
\alias{clDataCons,missing-method}
\title{Class "clDataCons"}
\description{This class is similar to \link{clDataVal} but stores the information aggregated at space,time, technical level. Normaly it should be post-stratified information that can be used to apply statistical methods safely (we hope).}
\section{Objects from the Class}{
The creator function "clDataCons" can be called to create objects from this class.
}
\section{Slots}{
	 \describe{
    \item{\code{desc}:}{Object of class \code{"character"} describing the object.}
    \item{\code{cl}:}{Object of class \code{"data.frame"} with effort information }
  }
}
\section{Methods}{
No methods defined with class "clDataCons" in the signature.
}
\author{Ernesto Jardim <ernesto@ipimar.pt>}
\examples{
showClass("clDataCons")
}
\keyword{classes}

