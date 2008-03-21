\name{DeltaID-class}
\docType{class}
\alias{DeltaID}
\alias{DeltaID-class}
\title{Class "DeltaID"}
\description{Any object of this class is created by \emph{plot.Delta} procedure.}

\section{Slots}{
\tabular{lrll}{
\bold{slot} \tab  \tab \bold{class} \tab \bold{description} \cr
\bold{\code{OutP}} \tab \tab \code{list} \tab Data from \code{plot.Delta} identification tool \cr
 \tab species \tab \code{character} \tab Specified species. \cr
 \tab tabID \tab \code{character} \tab Data from identified samples. \cr
 \tab tab \tab \code{list} \tab Resulting data from \code{Delta} procedure.  
}}

\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "DeltaID")}: graphical display of length distributions for identified samples.}
    \item{plot.Samp}{\code{signature(x = "DeltaID")}: graphical display of length distribution for one given sample, compared with overall relative length distribution (individualized \emph{plot} procedure).}
	 }
}


\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}
\author{Mathieu Merzereaud}
\seealso{\code{\link{Delta}}, \code{\link{DeltA-class}}, \code{\link{plot.Delta}}, \code{\link{plot.DeltaID}}, \code{\link{plot.LD}}, \code{\link{plot.Samp}}
}

\examples{
showClass("DeltaID")
}
\keyword{classes}

