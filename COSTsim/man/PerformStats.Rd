\name{PerformStats}
\alias{PerformStats}
\alias{PerformStats-class}
\docType{class}
\title{Class "PerformStats"}
\description{
  The PerformStats class stores the performance statistics to evaluate accuracy, bias and precision of different estimators. It is the outcome 
  object after applying the function \code{PerformStats}. 
}

\section{Slots}{
  \tabular{lrl}{
  \bold{slot} \tab \bold{class} \tab \bold{description} \cr
  \bold{\code{desc}} \tab \code{character} \tab object description \cr
  \bold{\code{species}} \tab \code{character} \tab species description (recall of SL\$spp (+ SL\$sex)) \cr
  \bold{\code{catchCat}} \tab \code{character} \tab recall of the catch category (discards/landings) \cr
  \bold{\code{param}} \tab \code{character} \tab recall of the parameter estimated (N, W, maturity, sex-ratio,...) \cr
	\bold{\code{strataDesc}} \tab \code{strIni} \tab time, space and technical stratification considered \cr
	\bold{\code{methodDesc}} \tab \code{character} \tab recall of the method (analytical, bootstrap, bayesian) \cr
	\bold{\code{nSamples}} \tab \code{numeric} \tab number of samples \cr
  \bold{\code{ageTrue}} \tab \code{data.frame} \tab True age structure (param-at-length) \cr
  \bold{\code{ageEst}} \tab \code{data.frame} \tab mean estimates of the age structure (param-at-length) \cr
	\bold{\code{ageAcc}} \tab \code{list} \tab a list with Accuracy statistics at age \cr
	\bold{\code{ageBias}} \tab \code{list} \tab a list with Bias statistics at age \cr
	\bold{\code{agePrec}} \tab \code{list} \tab a list with Precision statistics at age \cr
  \bold{\code{lenTrue}} \tab \code{data.frame} \tab True age structure (param-at-length) \cr
  \bold{\code{lenEst}} \tab \code{data.frame} \tab mean estimates of the length structure (param-at-length) \cr
	\bold{\code{lenAcc}} \tab \code{list} \tab a list with Accuracy statistics at length \cr
	\bold{\code{lenBias}} \tab \code{list} \tab a list with Bias statistics at length \cr
	\bold{\code{lenPrec}} \tab \code{list} \tab a list with Precision statistics at length \cr
  \bold{\code{totalNTrue}} \tab \code{data.frame} \tab True age structure (param-at-length) \cr
  \bold{\code{totalNEst}} \tab \code{data.frame} \tab mean estimates of the length structure (param-at-length) \cr
	\bold{\code{totalNAcc}} \tab \code{list} \tab a list with Accuracy statistics at length \cr
	\bold{\code{totalNBias}} \tab \code{list} \tab a list with Bias statistics at length \cr
	\bold{\code{totalNPrec}} \tab \code{list} \tab a list with Precision statistics at length \cr
	\bold{\code{totalWTrue}} \tab \code{data.frame} \tab True age structure (param-at-length) \cr
  \bold{\code{totalWEst}} \tab \code{data.frame} \tab mean estimates of the length structure (param-at-length) \cr
	\bold{\code{totalWAcc}} \tab \code{list} \tab a list with Accuracy statistics at length \cr
	\bold{\code{totalWBias}} \tab \code{list} \tab a list with Bias statistics at length \cr
	\bold{\code{totalWPrec}} \tab \code{list} \tab a list with Precision statistics at length \cr
  }
  
  The accuracy, bias and precision statistics are described in detail in the COST manual
  
}
\author{Dorleta Garcia \email{dgarcia@azti.es}}

\examples{
showClass("PerformStats")
}

\keyword{classes}
