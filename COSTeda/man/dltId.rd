\name{dltId-class}
\docType{class}
\alias{dltId}
\alias{dltId-class}
\title{Class "dltId"}
\description{Any object of this class is created by \emph{dltPlot} procedure.}

\section{Slots}{
\tabular{lrll}{
\bold{slot} \tab  \tab \bold{class} \tab \bold{description} \cr
\bold{\code{outP}} \tab \tab \code{list} \tab Data from \code{dltPlot} identification tool \cr
 \tab species \tab \code{character} \tab Specified species. \cr
 \tab sampId \tab \code{data.frame} \tab Specification of identified samples. \cr
 \tab tabId \tab \code{data.frame} \tab Length distribution data from identified samples. \cr
 \tab tab \tab \code{data.frame} \tab Resulting data from \code{dltCalc} procedure.  
}}

\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "dltId")}: graphical display of length distributions for identified samples, compared with overall relative length distribution.}
    \item{smpPlot}{\code{signature(x = "dltId")}: graphical display of length distribution for one given sample, compared with overall relative length distribution (individualized \emph{plot} procedure).}
	 }
}


\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}
\author{Mathieu Merzereaud}
\seealso{\code{\link{dltCalc}}, \code{\link{dltCls}}, \code{\link{dltPlot}}, \code{\link{plot.dltId}}, \code{\link{smpPlot}}, \code{\link{lenDisPlot}}
}

\examples{
showClass("dltId")
}
\keyword{classes}

