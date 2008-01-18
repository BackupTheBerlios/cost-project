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
\tabular{lrll}{
\bold{slot} \tab \bold{columns} \tab \bold{class} \tab \bold{description} \cr
\bold{\code{desc}} \tab \tab \code{character} \tab object description \cr
\bold{\code{cl}} \tab \tab \code{data.frame} \tab landings information \cr
 \tab landCtry \tab \code{factor} \tab landing country \cr
 \tab vslFlgCtry \tab \code{factor} \tab vessel flag country \cr
 \tab year \tab \code{numeric} \tab year \cr
 \tab quarter \tab \code{numeric} \tab quarter \cr 
 \tab month \tab \code{numeric} \tab month \cr
 \tab area \tab \code{factor} \tab area \cr
 \tab rect \tab \code{factor} \tab statistical rectangle \cr 
 \tab spp \tab \code{factor} \tab species \cr 
 \tab landCat \tab \code{factor} \tab landing category \cr 
 \tab commCatScl \tab \code{factor} \tab commercial category scale \cr
 \tab commCat \tab \code{factor} \tab commercial category \cr
 \tab foCatNat \tab \code{factor} \tab fishing operations category, national \cr
 \tab foCatEu5 \tab \code{factor} \tab fishing operations category, european level 5 \cr
 \tab foCatEu6 \tab \code{factor} \tab fishing operations category, european level 6 \cr
 \tab unallocCatchWt \tab \code{numeric} \tab unallocated catch weight \cr
 \tab misRepCatchWt \tab \code{numeric} \tab area misreported catch weight \cr
 \tab landWt \tab \code{numeric} \tab official landings weight \cr
 \tab landMult \tab \code{numeric} \tab landings multiplier \cr
 \tab landValue \tab \code{numeric} \tab official landings value \cr
}
}
\section{Methods}{
  \describe{
    \item{cl}{\code{signature(object = "clData")}: ... }
    \item{clDataVal}{\code{signature(object = "clData")}: ... }
    \item{desc}{\code{signature(object = "clData")}: ... }
    \item{dim}{\code{signature(x = "clData")}: ... }
    \item{head}{\code{signature(x = "clData")}: ... }
    \item{rbind2}{\code{signature(x = "clData", y = "clData")}: ... }
    \item{subset}{\code{signature(x = "clData")}: ... }
    \item{summary}{\code{signature(object = "clData")}: ... }
    \item{tail}{\code{signature(x = "clData")}: ... }
	 }
}
\author{Ernesto Jardim <ernesto@ipimar.pt>}
\examples{
showClass("clData")
}
\keyword{classes}
