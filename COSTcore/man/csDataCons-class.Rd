\name{csDataCons-class}
\docType{class}
\alias{csDataCons}
\alias{csDataCons-class}
\alias{csDataCons-methods}
\alias{csDataCons,csDataVal-method}
\alias{csDataCons,missing-method}
\title{Class "csDataCons"}
\description{This class is similar to \link{csDataVal} but stores the information aggregated at space,time, technical level. Normaly it should be post-stratified information that can be used to apply statistical methods safely (we hope).}
\section{Objects from the Class}{
The creator function "csDataCons" can be called to create objects from this class.
}
\section{Slots}{
\tabular{lrll}{
\bold{slot} \tab \bold{columns} \tab \bold{class} \tab \bold{description} \cr
\bold{\code{desc}} \tab \tab \code{character} \tab object description \cr
 \tab \tab \tab \cr
\bold{\code{tr}} \tab \tab \code{data.frame} \tab  trip information \cr
 \tab PSUid \tab \code{numeric} \tab primary sampling unit id \cr
 \tab sampType \tab \code{factor} \tab sampling type \cr
 \tab landCtry \tab \code{factor} \tab landing country \cr
 \tab vslFlgCtry \tab \code{factor} \tab vessel flag country \cr
 \tab year \tab \code{numeric} \tab year \cr
 \tab proj \tab \code{factor} \tab project \cr
 \tab trpCode \tab \code{factor} \tab trip code \cr
 \tab vslLen \tab \code{numeric} \tab vessel length-over-all \cr 
 \tab vslPwr \tab \code{numeric} \tab vessel power (kw)\cr 
 \tab vslSize \tab \code{numeric} \tab vessel size (GRT)\cr 
 \tab vsType \tab \code{factor} \tab vessel type \cr 
 \tab foNum \tab \code{numeric} \tab number of fishing operations \cr 
 \tab daysAtSea \tab \code{numeric} \tab effort in days at sea \cr 
 \tab vslId \tab \code{numeric} \tab vessel identifier \cr 
 \tab sampCtry \tab \code{factor} \tab sampling country \cr 
 \tab sampMeth \tab \code{factor} \tab sampling method \cr
 \tab \tab \tab \cr
\bold{\code{hh}} \tab \tab \code{data.frame} \tab haul information \cr
 \tab PSUid \tab \code{numeric} \tab primary sampling unit id \cr
 \tab staNum \tab \code{numeric} \tab station number \cr
 \tab SSUid \tab \code{numeric} \tab secondary sampling unit id \cr
 \tab time \tab \code{factor} \tab fishing operation validity \cr
 \tab space \tab \code{factor} \tab date \cr
 \tab technical \tab \code{factor} \tab date \cr
 \tab aggLev \tab \code{factor} \tab aggregation level \cr
 \tab foDur \tab \code{numeric} \tab fishing operation duration \cr
 \tab \tab \tab \cr
\bold{\code{sl}} \tab \tab \code{data.frame} \tab species list \cr
 \tab PSUid \tab \code{numeric} \tab primary sampling unit id \cr
 \tab SSUid \tab \code{numeric} \tab secondary sampling unit id \cr
 \tab spp \tab \code{factor} \tab species \cr 
 \tab catchCat \tab \code{factor} \tab catch category \cr 
 \tab landCat \tab \code{factor} \tab landing category \cr 
 \tab commCatScl \tab \code{factor} \tab commercial category scale \cr
 \tab commCat \tab \code{factor} \tab commercial category \cr
 \tab subSampCat \tab \code{factor} \tab subsampling category \cr
 \tab TSUid \tab \code{numeric} \tab terciary sampling unit id \cr
 \tab time \tab \code{factor} \tab fishing operation validity \cr
 \tab space \tab \code{factor} \tab date \cr
 \tab technical \tab \code{factor} \tab date \cr
 \tab sort \tab \code{factor} \tab validity code \cr 
 \tab wt \tab \code{numeric} \tab weight \cr 
 \tab subSampWt \tab \code{numeric} \tab subsample weight \cr 
 \tab lenCode \tab \code{factor} \tab length code \cr
 \tab \tab \tab \cr
\bold{\code{hl}} \tab \tab \code{data.frame} \tab length frequencies \cr
 \tab PSUid \tab \code{numeric} \tab primary sampling unit id \cr
 \tab SSUid \tab \code{numeric} \tab secondary sampling unit id \cr
 \tab TSUid \tab \code{numeric} \tab terciary sampling unit id \cr
 \tab time \tab \code{factor} \tab fishing operation validity \cr
 \tab space \tab \code{factor} \tab date \cr
 \tab technical \tab \code{factor} \tab date \cr
 \tab sort \tab \code{factor} \tab sampling type \cr
 \tab lenCls \tab \code{numeric} \tab length class \cr
 \tab lenNum \tab \code{numeric} \tab frequency at length \cr
 \tab \tab \tab \cr
\bold{\code{ca}} \tab \tab \code{data.frame} \tab growth and maturity \cr
 \tab PSUid \tab \code{numeric} \tab primary sampling unit id \cr
 \tab SSUid \tab \code{numeric} \tab secondary sampling unit id \cr
 \tab time \tab \code{factor} \tab fishing operation validity \cr
 \tab space \tab \code{factor} \tab date \cr
 \tab technical \tab \code{factor} \tab date \cr
 \tab sort \tab \code{factor} \tab sampling type \cr
 \tab spp \tab \code{factor} \tab species \cr 
 \tab sex \tab \code{factor} \tab gender \cr 
 \tab catchCat \tab \code{factor} \tab catch category \cr 
 \tab landCat \tab \code{factor} \tab landing category \cr 
 \tab commCatScl \tab \code{factor} \tab commercial category scale \cr
 \tab commCat \tab \code{factor} \tab commercial category \cr
 \tab stock \tab \code{factor} \tab stock \cr
 \tab area \tab \code{factor} \tab area \cr
 \tab rect \tab \code{factor} \tab statistical rectangle \cr
 \tab lenCls \tab \code{numeric} \tab length class \cr
 \tab age \tab \code{numeric} \tab age \cr
 \tab fishId \tab \code{numeric} \tab individual identification \cr
 \tab lenCode \tab \code{factor} \tab length code \cr
 \tab plusGrp \tab \code{factor} \tab plus group \cr
 \tab otoWt \tab \code{numeric} \tab otolith weight \cr
 \tab otoSide \tab \code{factor} \tab otolith side \cr
 \tab indWt \tab \code{numeric} \tab individual weight \cr
 \tab matScale \tab \code{factor} \tab maturity scale \cr
 \tab matStage \tab \code{factor} \tab maturity stage \cr
}
}
\section{Methods}{
  \describe{
    \item{ca}{\code{signature(object = "csDataCons")}: ... }
    \item{desc}{\code{signature(object = "csDataCons")}: ... }
    \item{dim}{\code{signature(x = "csDataCons")}: ... }
    \item{head}{\code{signature(x = "csDataCons")}: ... }
    \item{hh}{\code{signature(object = "csDataCons")}: ... }
    \item{hl}{\code{signature(object = "csDataCons")}: ... }
    \item{rbind2}{\code{signature(x = "csDataCons", y = "csDataCons")}: ... }
    \item{sl}{\code{signature(object = "csDataCons")}: ... }
    \item{subset}{\code{signature(x = "csDataCons")}: ... }
    \item{summary}{\code{signature(object = "csDataCons")}: ... }
    \item{tail}{\code{signature(x = "csDataCons")}: ... }
    \item{tr}{\code{signature(object = "csDataCons")}: ... }
	 }
}
\author{Ernesto Jardim <ernesto@ipimar.pt>}
\examples{
showClass("csDataCons")
}
\keyword{classes}

