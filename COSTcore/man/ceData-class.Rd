\name{ceData-class}
\docType{class}
\alias{ceData}
\alias{ceData-class}
\alias{ceData-methods}
\alias{ceData,data.frame-method}
\alias{ceData,matrix-method}
\alias{ceData,character-method}
\alias{ceData,missing-method}
\title{Class "ceData"}
\description{This class stores commercial fisheries effort information and maps the CE tables in Fish Frame version 5.}
\section{Objects from the Class}{
The creator function "ceData" can be called to create objects from this class.
}
\section{Slots}{
\tabular{lrll}{
\bold{slot} \tab \tab \bold{class} \tab \bold{description} \cr
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
No methods defined with class "csData" in the signature.
}
\author{Ernesto Jardim <ernesto@ipimar.pt>}
\examples{
showClass("ceData")
}
\keyword{classes}
