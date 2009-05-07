\name{trueData}
\alias{trueData}
\docType{class}
\title{Class "trueData"}
\description{
  The trueData class stores the true data set
}

\section{Objects from the Class}{
The creator function \code{trueData} can be called to create objects from this class.
}

\section{Slots}{
  \tabular{lrl}{
  \bold{slot} \tab \bold{class} \tab \bold{description} \cr
  \bold{\code{desc}} \tab \code{character} \tab object description \cr
  \bold{\code{species}} \tab \code{character} \tab species description (recall of SL\$spp (+ SL\$sex)) \cr
  \bold{\code{strataDesc}} \tab \code{StrIni} \tab time, space and technical stratification considered \cr
  \bold{\code{lal}} \tab \code{data.frame} \tab true landings at age \cr
  \bold{\code{laa}} \tab \code{data.frame} \tab true landings at length \cr
  \bold{\code{dal}} \tab \code{data.frame} \tab true discards at length \cr
  \bold{\code{dtw}} \tab \code{data.frame} \tab true discards total weight \cr  
  }
}

\author{Dorleta García \email{dgarcia@azti.es}}

\examples{
showclass("trueData")
}

\keyword{classes}
