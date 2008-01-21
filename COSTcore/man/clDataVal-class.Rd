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
\tabular{lrll}{
\bold{slot} \tab \bold{columns} \tab \bold{class} \tab \bold{description} \cr
\bold{\code{desc}} \tab \tab \code{character} \tab object description \cr
\bold{\code{cl}} \tab \tab \code{data.frame} \tab landings information \cr
 \tab landCtry \tab \code{character} \tab landing country \cr
 \tab vslFlgCtry \tab \code{character} \tab vessel flag country \cr
 \tab year \tab \code{numeric} \tab year \cr
 \tab quarter \tab \code{numeric} \tab quarter \cr 
 \tab month \tab \code{numeric} \tab month \cr
 \tab area \tab \code{character} \tab area \cr
 \tab rect \tab \code{character} \tab statistical rectangle \cr 
 \tab spp \tab \code{character} \tab species \cr 
 \tab landCat \tab \code{character} \tab landing category \cr 
 \tab commCatScl \tab \code{character} \tab commercial category scale \cr
 \tab commCat \tab \code{character} \tab commercial category \cr
 \tab foCatNat \tab \code{character} \tab fishing operations category, national \cr
 \tab foCatEu5 \tab \code{character} \tab fishing operations category, european level 5 \cr
 \tab foCatEu6 \tab \code{character} \tab fishing operations category, european level 6 \cr
 \tab unallocCatchWt \tab \code{numeric} \tab unallocated catch weight \cr
 \tab misRepCatchWt \tab \code{numeric} \tab area misreported catch weight \cr
 \tab landWt \tab \code{numeric} \tab official landings weight \cr
 \tab landMult \tab \code{numeric} \tab landings multiplier \cr
 \tab landValue \tab \code{numeric} \tab official landings value \cr
}
}
\section{Methods}{

  \describe{
    \item{cl}{\code{signature(object = "clData")}: \code{cl} slot accessor.}
    \item{clDataCons}{\code{signature(object = "clDataVal")}: creator of Consolidated data class.}
    \item{desc}{\code{signature(object = "clDataVal")}: \code{desc} slot accessor.}
    \item{dim}{\code{signature(x = "clDataVal")}: dimensions of \code{cl} slot (see \code{\link{dim}}).}
    \item{head}{\code{signature(x = "clDataVal")}: first rows of the \code{cl} slot (see \code{\link{head}}).}
    \item{rbind2}{\code{signature(x = "clDataVal", y = "clDataVal")}: bind two \code{clDataVal} objects by adding rows (see \code{\link{rbind2}}).}
    \item{subset}{\code{signature(x = "clDataVal")}: subset a \code{clDataVal} object (see \code{\link{subset}}).}
    \item{summary}{\code{signature(object = "clDataVal")}: summary a \code{clDataVal} object (see \code{\link{summary}}).}
    \item{tail}{\code{signature(x = "clDataVal")}: last rows of the \code{cl} slot (see \code{\link{tail}}).}
	 }
}
\author{Ernesto Jardim <ernesto@ipimar.pt>}
\examples{
showClass("clDataVal")
}
\keyword{classes}

