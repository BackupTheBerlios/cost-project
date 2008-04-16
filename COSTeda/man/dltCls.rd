\name{dltCls-class}
\docType{class}
\alias{dltCls}
\alias{dltCls-class}
\title{Class "dltCls"}
\description{Any object of this class is created by \emph{dltCalc} procedure.}

\section{Slots}{
\tabular{lrll}{
\bold{slot} \tab  \tab \bold{class} \tab \bold{description} \cr
\bold{\code{strPar}} \tab \tab \code{list} \tab Informations from \code{dltCalc} function call \cr
 \tab species \tab \code{character} \tab Specified species. \cr
 \tab timeStrata \tab \code{character} \tab Specified time stratification. \cr
 \tab spaceStrata \tab \code{character} \tab Specified space stratification. \cr
 \tab techStrata \tab \code{character} \tab Specified technical stratification. \cr
\bold{\code{outP}} \tab \tab \code{list} \tab Numerical outputs from \code{Delta} method \cr
 \tab DeltaMatrix \tab \code{character} \tab Sum of squared delta values of each sample within  \cr
                              \tab \tab \tab each strata and length class (only if \code{indSamp=FALSE} \cr
                              \tab \tab \tab in \code{dltCalc} function call).\cr
 \tab NkMatrix \tab \code{character} \tab Number of samples within each strata (only if \code{indSamp=FALSE} \cr 
                           \tab \tab \tab in \code{dltCalc} function call). \cr
 \tab WkMatrix \tab \code{character} \tab Sampled weight within each strata in grams (only if \code{indSamp=FALSE} \cr 
                           \tab \tab \tab in \code{dltCalc} function call). \cr
 \tab SampDeltaMat \tab \code{character} \tab Data.frame of delta values within each sample \cr 
                               \tab \tab \tab (only if \code{indSamp=TRUE} in \code{dltCalc} function call). \cr
 \tab tab \tab \code{character} \tab Data.frame resulting from treatment (only if \cr 
                      \tab \tab \tab \code{indSamp=TRUE} in \code{dltCalc} function call). \cr
 \tab DFsamp \tab \code{character} \tab Informations about each sample (only if \code{indSamp=TRUE} \cr
                      \tab \tab \tab in \code{dltCalc} function call).
}
}


\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}
\author{Mathieu Merzereaud}
\seealso{\code{\link{dltCalc}}, \code{\link{dltId}}, \code{\link{dltPlot}}, \code{\link{plot.dltId}}, \code{\link{smpPlot}}, \code{\link{lenDisPlot}}
}

\examples{
showClass("dltCls")
}
\keyword{classes}

