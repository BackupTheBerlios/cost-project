\name{clDataVal-class}
\docType{class}
\alias{clDataVal}
\alias{clDataVal-class}
\alias{clDataVal-methods}
\alias{clDataVal,clData-method}
\alias{clDataVal,missing-method}
\title{Class "clDataVal"}
\description{This class is equal to \link{clData} but stores the information after being processed for coding errors, typos, outliers, etc. The information contained is considered to be validated, which means it can now be used for post-stratifications and other aggregation exercises.}
\section{Objects from the Class}{
The creator function "clDataVal" can be called to create objects from this class.
}
\section{Slots}{
	 \describe{
    \item{\code{desc}:}{Object of class \code{"character"} describing the object.}
    \item{\code{cl}:}{Object of class \code{"data.frame"} with effort information }
  }
}
\section{Methods}{
No methods defined with class "clDataVal" in the signature.
}
\author{Ernesto Jardim <ernesto@ipimar.pt>}
\examples{
showClass("clDataVal")
}
\keyword{classes}

