\name{csDataVal-class}
\docType{class}
\alias{csDataVal}
\alias{csDataVal-class}
\alias{csDataVal-methods}
\alias{csDataVal,csData-method}
\alias{csDataVal,missing-method}
\title{Class "csDataVal"}
\description{This class is equal to \link{csData} but stores the information after being processed for coding errors, typos, outliers, etc. The information contained is considered to be validated, which means it can now be used for post-stratifications and other aggregation exercises.}
\section{Objects from the Class}{
The creator function "csDataVal" can be called to create objects from this class.
}
\section{Slots}{
	 \describe{
    \item{\code{desc}:}{Object of class \code{"character"} describing the object.}
    \item{\code{tr}:}{Object of class \code{"data.frame"} with trip information.}
    \item{\code{hh}:}{Object of class \code{"data.frame"} with haul information.}
    \item{\code{sl}:}{Object of class \code{"data.frame"} with species list information.}
    \item{\code{hl}:}{Object of class \code{"data.frame"} with length frequencies information.}
    \item{\code{ca}:}{Object of class \code{"data.frame"} with growth and maturity information.}
  }
}
\section{Methods}{
No methods defined with class "csDataVal" in the signature.
}
\author{Ernesto Jardim <ernesto@ipimar.pt>}
\examples{
showClass("csDataVal")
}
\keyword{classes}
