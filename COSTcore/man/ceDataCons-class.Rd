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
\tabular{lrll}{
\bold{slot} \tab \bold{columns} \tab \bold{class} \tab \bold{description} \cr
\bold{\code{desc}} \tab \tab \code{character} \tab object description \cr
\bold{\code{ce}} \tab \tab \code{data.frame} \tab effort information \cr
 \tab vslFlgCtry \tab \code{factor} \tab vessel flag country \cr
 \tab time \tab \code{factor} \tab year \cr
 \tab space \tab \code{factor} \tab quarter \cr 
 \tab technical \tab \code{factor} \tab month \cr
 \tab trpNum \tab \code{numeric} \tab number of trips \cr
 \tab foNum \tab \code{numeric} \tab number of fishing operations \cr
 \tab foDur \tab \code{numeric} \tab duration of fishing operations \cr
 \tab effKwDays \tab \code{numeric} \tab effort in kw days \cr
 \tab effGtDays \tab \code{numeric} \tab effort in T days \cr
 \tab daysAtSea \tab \code{numeric} \tab effort in days at sea \cr	
}
}
\section{Methods}{
  \describe{
    \item{ce}{\code{signature(object = "ceDataCons")}: ... }
    \item{desc}{\code{signature(object = "ceDataCons")}: ... }
    \item{dim}{\code{signature(x = "ceDataCons")}: ... }
    \item{head}{\code{signature(x = "ceDataCons")}: ... }
    \item{rbind2}{\code{signature(x = "ceDataCons", y = "ceDataCons")}: ... }
    \item{subset}{\code{signature(x = "ceDataCons")}: ... }
    \item{summary}{\code{signature(object = "ceDataCons")}: ... }
    \item{tail}{\code{signature(x = "ceDataCons")}: ... }
	 }
}
\author{Ernesto Jardim <ernesto@ipimar.pt>}
\examples{
showClass("ceDataCons")
}
\keyword{classes}

