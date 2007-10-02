\name{ceData-class}
\docType{class}
\alias{ceData}
\alias{ceData-class}
\alias{ceData-methods}
\alias{ceData,data.frame-method}
\alias{ceData,character-method}
\alias{ceData,missing-method}
\title{Class "ceData"}
\description{This class stores commercial fisheries effort information and maps the CE tables in Fish Frame version 5.}
\section{Objects from the Class}{
The creator function "ceData" can be called to create objects from this class.
}
\section{Slots}{
	 \describe{
    \item{\code{desc}:}{Object of class \code{"character"} describing the object.}
    \item{\code{ce}:}{Object of class \code{"data.frame"} with effort information }
  }
}
\section{Methods}{
No methods defined with class "csData" in the signature.
}
\author{Ernesto Jardim <ernesto@ipimar.pt>}
\examples{
showClass("ceData")
}
\keyword{classes}
