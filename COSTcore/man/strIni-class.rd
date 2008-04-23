\name{strIni-class}
\docType{class}
\alias{strIni}
\alias{strIni-class}
\alias{strIni-methods}
\title{Class "strIni"}
\description{The strIni class stores stratification and recoding information that is required for consolidated objects creation process.}
\section{Objects from the Class}{
The creator function "strIni" can be called to create objects from this class. Default values for slots are NAs, that means no stratification or no recoding)
}
\section{Slots}{
\tabular{lrll}{
\bold{slot} \tab \bold{columns} \tab \bold{class} \tab \bold{description} \cr
\bold{\code{timeStrata}} \tab \tab \code{character} \tab time stratification (e.g "quarter", "month",...)\cr
\bold{\code{spaceStrata}} \tab \tab \code{character} \tab space stratification (e.g "area",...)\cr
\bold{\code{techStrata}} \tab \tab \code{character} \tab technical stratification (e.g "commCat", "foCatEu5",...)\cr
\bold{\code{tpRec}} \tab \tab \code{list} \tab list for time strata recoding \cr
\tab \tab \tab (e.g \code{list(from="1",to="2")}) \cr
\bold{\code{spRec}} \tab \tab \code{list} \tab list for space strata recoding \cr
\tab \tab \tab (e.g \code{list(from=c("7D","7D1"),to=c("27.7.d","27.7.d"))}) \cr
\bold{\code{tcRec}} \tab \tab \code{list} \tab list for technical strata recoding \cr
\tab \tab \tab (e.g \code{list(from=c("OTB-DEF","OTB-MOL"),to=c("OTB","OTB"))}) 
}
}

\section{Methods}{

  \describe{
    \item{csDataCons}{\code{signature(object="csDataVal",objStrat="strIni")}:\code{csDataCons} class creator.}
    \item{clDataCons}{\code{signature(object="clDataVal",objStrat="strIni")}:\code{clDataCons} class creator.}
    \item{ceDataCons}{\code{signature(object="ceDataVal",objStrat="strIni")}:\code{ceDataCons} class creator.}
	 }
}
\author{Mathieu Merzereaud <Mathieu.Merzereaud@ifremer.fr>}
\examples{
showClass("strIni")
}
\keyword{classes}

