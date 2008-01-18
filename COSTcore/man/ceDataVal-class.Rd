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
\tabular{lrll}{
\bold{slot} \tab \bold{columns} \tab \bold{class} \tab \bold{description} \cr
\bold{\code{desc}} \tab \tab \code{character} \tab object description \cr
\bold{\code{ce}} \tab \tab \code{data.frame} \tab effort information \cr
 \tab vslFlgCtry \tab \code{factor} \tab vessel flag country \cr
 \tab year \tab \code{numeric} \tab year \cr
 \tab quarter \tab \code{numeric} \tab quarter \cr 
 \tab month \tab \code{numeric} \tab month \cr
 \tab area \tab \code{factor} \tab area \cr
 \tab rect \tab \code{factor} \tab statistical rectangle \cr 
 \tab foCatNat \tab \code{factor} \tab fishing operations category, national \cr
 \tab foCatEu5 \tab \code{factor} \tab fishing operations category, european level 5 \cr
 \tab foCatEu6 \tab \code{factor} \tab fishing operations category, european level 6 \cr
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
    \item{ce}{\code{signature(object = "ceData")}: ... }
    \item{ceDataVal}{\code{signature(object = "ceData")}: ... }
    \item{desc}{\code{signature(object = "ceData")}: ... }
    \item{dim}{\code{signature(x = "ceData")}: ... }
    \item{head}{\code{signature(x = "ceData")}: ... }
    \item{rbind2}{\code{signature(x = "ceData", y = "ceData")}: ... }
    \item{subset}{\code{signature(x = "ceData")}: ... }
    \item{summary}{\code{signature(object = "ceData")}: ... }
    \item{tail}{\code{signature(x = "ceData")}: ... }
	 }
}
\author{Ernesto Jardim <ernesto@ipimar.pt>}
\examples{
showClass("ceDataVal")
}
\keyword{classes}

