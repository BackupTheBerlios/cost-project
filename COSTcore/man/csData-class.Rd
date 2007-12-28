\name{csData-class}
\docType{class}
\alias{csData}
\alias{csData-class}
\alias{csData-methods}
\alias{csData,missing,missing,missing,missing,missing-method}
\alias{csData,data.frame,missing,missing,missing,missing-method}
\alias{csData,data.frame,data.frame,missing,missing,missing-method}
\alias{csData,data.frame,data.frame,data.frame,missing,missing-method}
\alias{csData,data.frame,data.frame,data.frame,data.frame,missing-method}
\alias{csData,data.frame,data.frame,data.frame,data.frame,data.frame-method}
\alias{csData,character,character,character,character,missing-method}
\alias{csData,character,character,character,character,character-method}
\title{Class "csData"}
\description{This class stores commercial fisheries sampling information and maps the CS tables in Fish Frame version 5.}
\section{Objects from the Class}{
The creator function "csData" can be called to create objects from this class.
}
\section{Slots}{
\tabular{lrll}{
\bold{slot} \tab \tab \bold{class} \tab \bold{description} \cr
\bold{\code{desc}} \tab \tab \code{character} \tab object description \cr
 \tab \tab \tab \cr
\bold{\code{tr}} \tab \tab \code{data.frame} \tab  trip information \cr
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
 \tab sampType \tab \code{factor} \tab sampling type \cr
 \tab landCtry \tab \code{factor} \tab landing country \cr
 \tab vslFlgCtry \tab \code{factor} \tab vessel flag country \cr
 \tab year \tab \code{numeric} \tab year \cr
 \tab proj \tab \code{factor} \tab project \cr
 \tab trpCode \tab \code{factor} \tab trip code \cr
 \tab staNum \tab \code{numeric} \tab station number \cr
 \tab foVal \tab \code{factor} \tab fishing operation validity \cr
 \tab aggLev \tab \code{factor} \tab aggregation level \cr
 \tab date \tab \code{factor} \tab date \cr
 \tab time \tab \code{factor} \tab time \cr
 \tab foDur \tab \code{numeric} \tab fishing operation duration \cr
 \tab latIni \tab \code{numeric} \tab latitude of shooting position in decimal degrees \cr
 \tab lonIni \tab \code{numeric} \tab longitude of shooting position in decimal degrees \cr
 \tab latFin \tab \code{numeric} \tab latitude of hauling position in decimal degrees \cr
 \tab lonFin \tab \code{numeric} \tab longitude of hauling position in decimal degrees \cr
 \tab area \tab \code{factor} \tab area \cr
 \tab rect \tab \code{factor} \tab statistical rectangle \cr
 \tab foDep \tab \code{numeric} \tab fishing operation depth \cr
 \tab waterDep \tab \code{numeric} \tab water depth \cr
 \tab foCatNat \tab \code{factor} \tab fishing operations category, national \cr
 \tab foCatEu5 \tab \code{factor} \tab fishing operations category, european level 5 \cr
 \tab foCatEu6 \tab \code{factor} \tab fishing operations category, european level 6 \cr
 \tab gear \tab \code{factor} \tab gear \cr
 \tab meshSize \tab \code{numeric} \tab mesh size \cr
 \tab selDev \tab \code{factor} \tab selectivity device \cr
 \tab meshSizeSelDev \tab \code{numeric} \tab mesh size in selectivity device \cr
 \tab \tab \tab \cr
\bold{\code{sl}} \tab \tab \code{data.frame} \tab species list \cr
 \tab sampType \tab \code{factor} \tab sampling type \cr
 \tab landCtry \tab \code{factor} \tab landing country \cr
 \tab vslFlgCtry \tab \code{factor} \tab vessel flag country \cr
 \tab year \tab \code{numeric} \tab year \cr
 \tab proj \tab \code{factor} \tab project \cr
 \tab trpCode \tab \code{factor} \tab trip code \cr
 \tab staNum \tab \code{numeric} \tab station number \cr
 \tab spp \tab \code{factor} \tab species \cr 
 \tab catchCat \tab \code{factor} \tab catch category \cr 
 \tab landCat \tab \code{factor} \tab landing category \cr 
 \tab commCatScl \tab \code{factor} \tab commercial category scale \cr
 \tab commCat \tab \code{factor} \tab commercial category \cr
 \tab subSampCat \tab \code{factor} \tab subsampling category \cr
 \tab valCode \tab \code{factor} \tab validity code \cr 
 \tab wt \tab \code{numeric} \tab weight \cr 
 \tab subSampWt \tab \code{numeric} \tab subsample weight \cr 
 \tab lenCode \tab \code{factor} \tab length code \cr
 \tab \tab \tab \cr
\bold{\code{hl}} \tab \tab \code{data.frame} \tab length frequencies \cr
 \tab sampType \tab \code{factor} \tab sampling type \cr
 \tab landCtry \tab \code{factor} \tab landing country \cr
 \tab vslFlgCtry \tab \code{factor} \tab vessel flag country \cr
 \tab year \tab \code{numeric} \tab year \cr
 \tab proj \tab \code{factor} \tab project \cr
 \tab trpCode \tab \code{factor} \tab trip code \cr
 \tab staNum \tab \code{numeric} \tab station number \cr
 \tab spp \tab \code{factor} \tab species \cr 
 \tab sex \tab \code{factor} \tab gender \cr 
 \tab catchCat \tab \code{factor} \tab catch category \cr 
 \tab landCat \tab \code{factor} \tab landing category \cr 
 \tab commCatScl \tab \code{factor} \tab commercial category scale \cr
 \tab commCat \tab \code{factor} \tab commercial category \cr
 \tab subSampCat \tab \code{factor} \tab subsampling category \cr
 \tab lenCls \tab \code{numeric} \tab length class \cr
 \tab lenNum \tab \code{numeric} \tab frequency at length \cr
 \tab \tab \tab \cr
\bold{\code{ca}} \tab \tab \code{data.frame} \tab growth and maturity \cr
 \tab sampType \tab \code{factor} \tab sampling type \cr
 \tab landCtry \tab \code{factor} \tab landing country \cr
 \tab vslFlgCtry \tab \code{factor} \tab vessel flag country \cr
 \tab year \tab \code{numeric} \tab year \cr
 \tab quarter \tab \code{numeric} \tab quarter \cr
 \tab month \tab \code{numeric} \tab month \cr
 \tab proj \tab \code{factor} \tab project \cr
 \tab trpCode \tab \code{factor} \tab trip code \cr
 \tab staNum \tab \code{numeric} \tab station number \cr
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
No methods defined with class "csData" in the signature.
}
\author{Ernesto Jardim <ernesto@ipimar.pt>}
\examples{
showClass("csData")
}
\keyword{classes}
