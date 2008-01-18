\name{clDataCons-class}
\docType{class}
\alias{clDataCons}
\alias{clDataCons-class}
\alias{clDataCons-methods}
\alias{clDataCons,clDataVal-method}
\alias{clDataCons,missing-method}
\title{Class "clDataCons"}
\description{This class is similar to \link{clDataVal} but stores the information aggregated at space,time, technical level. Normaly it should be post-stratified information that can be used to apply statistical methods safely (we hope).}
\section{Objects from the Class}{
The creator function "clDataCons" can be called to create objects from this class.
}
\section{Slots}{
\tabular{lrll}{
\bold{slot} \tab \bold{columns} \tab \bold{class} \tab \bold{description} \cr
\bold{\code{desc}} \tab \tab \code{character} \tab object description \cr
\bold{\code{cl}} \tab \tab \code{data.frame} \tab landings information \cr
 \tab landCtry \tab \code{factor} \tab landing country \cr
 \tab vslFlgCtry \tab \code{factor} \tab vessel flag country \cr
 \tab time \tab \code{numeric} \tab year \cr
 \tab space \tab \code{numeric} \tab quarter \cr 
 \tab technocal \tab \code{numeric} \tab month \cr
 \tab spp \tab \code{factor} \tab species \cr 
 \tab landCat \tab \code{factor} \tab landing category \cr 
 \tab commCatScl \tab \code{factor} \tab commercial category scale \cr
 \tab commCat \tab \code{factor} \tab commercial category \cr
 \tab unallocCatchWt \tab \code{numeric} \tab unallocated catch weight \cr
 \tab misRepCatchWt \tab \code{numeric} \tab area misreported catch weight \cr
 \tab landWt \tab \code{numeric} \tab official landings weight \cr
 \tab landMult \tab \code{numeric} \tab landings multiplier \cr
 \tab landValue \tab \code{numeric} \tab official landings value \cr
}
}
\section{Methods}{
  \describe{
    \item{dim}{\code{signature(x = "clDataCons")}: ... }
    \item{head}{\code{signature(x = "clDataCons")}: ... }
    \item{rbind2}{\code{signature(x = "clDataCons", y = "clDataCons")}: ... }
    \item{subset}{\code{signature(x = "clDataCons")}: ... }
    \item{summary}{\code{signature(object = "clDataCons")}: ... }
    \item{tail}{\code{signature(x = "clDataCons")}: ... }
	 }
}
\author{Ernesto Jardim <ernesto@ipimar.pt>}
\examples{
showClass("clDataCons")
}
\keyword{classes}

