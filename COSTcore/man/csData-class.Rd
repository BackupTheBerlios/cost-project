\name{csData-class}
\docType{class}
\alias{csData}
\alias{csData-class}
\alias{csData-methods}
\alias{csData,data.frame,data.frame,data.frame,data.frame,missing-method}
\alias{csData,data.frame,data.frame,data.frame,data.frame,data.frame-method}
\alias{csData,character,character,character,character,missing-method}
\alias{csData,missing,missing,missing,missing,missing-method}
\title{Class "csData"}
\description{This class stores commercial fisheries sampling information and maps the CS tables in Fish Frame version 5.}
\section{Objects from the Class}{
The creator function "csData" can be called to create objects from this class.
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
No methods defined with class "csData" in the signature.
}
\author{Ernesto Jardim <ernesto@ipimar.pt>}
\examples{
showClass("csData")
}
\keyword{classes}
