\name{clData-class}
\docType{class}
\alias{clData}
\alias{clData-class}
\alias{clData-methods}
\alias{clData,data.frame-method}
\alias{clData,matrix-method}
\alias{clData,character-method}
\alias{clData,missing-method}
\title{Class "clData"}
\description{This class stores commercial fisheries landings information and maps the CL tables in Fish Frame version 5.}
\section{Objects from the Class}{
The creator function "clData" can be called to create objects from this class.
}
\section{Slots}{
	 \describe{
    \item{\code{desc}:}{Object of class \code{"character"} describing the object.}
    \item{\code{cl}:}{Object of class \code{"data.frame"} with landings information }
  }
}
\section{Methods}{
No methods defined with class "clData" in the signature.
}
\author{Ernesto Jardim <ernesto@ipimar.pt>}
\examples{
showClass("clData")
}
\keyword{classes}
