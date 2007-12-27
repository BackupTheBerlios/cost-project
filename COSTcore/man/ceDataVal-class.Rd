\name{ceDataVal-class}
\docType{class}
\alias{ceDataVal}
\alias{ceDataVal-class}
\alias{ceDataVal-methods}
\alias{ceDataVal,ceData-method}
\alias{ceDataVal,missing-method}
\title{Class "ceDataVal"}
\description{This class is equal to \link{ceData} but stores the information after being processed for coding errors, typos, outliers, etc. The information contained is considered to be validated, which means it can now be used for post-stratifications and other aggregation exercises.}
\section{Objects from the Class}{
The creator function "ceDataVal" can be called to create objects from this class.
}
\section{Slots}{
	 \describe{
    \item{\code{desc}:}{Object of class \code{"character"} describing the object.}
    \item{\code{ce}:}{Object of class \code{"data.frame"} with effort information }
  }
}
\section{Methods}{
No methods defined with class "ceDataVal" in the signature.
}
\author{Ernesto Jardim <ernesto@ipimar.pt>}
\examples{
showClass("ceDataVal")
}
\keyword{classes}

