\name{DeltA-class}
\docType{class}
\alias{DeltA}
\alias{DeltA-class}
\title{Class "DeltA"}
\description{Any object of this class is created by \emph{Delta} procedure.}

\section{Slots}{
\tabular{lrll}{
\bold{slot} \tab  \tab \bold{class} \tab \bold{description} \cr
\bold{\code{Parameters}} \tab \tab \code{list} \tab Informations about Delta function call \cr
 \tab species \tab \code{character} \tab Specified species (from \emph{Delta.list} object). \cr
 \tab tempStrata \tab \code{character} \tab Specified time stratification. \cr
 \tab spaceStrata \tab \code{character} \tab Specified space stratification. \cr
 \tab techStrata \tab \code{character} \tab Specified technical stratification. \cr
 \tab elmts \tab \code{character} \tab List of components within stratification. \cr
\bold{\code{OutP}} \tab \tab \code{list} \tab Numerical outputs from \code{Delta} method \cr
 \tab DeltaMatrix \tab \code{character} \tab Sum of squared delta values of each sample within each strata and length class (only if \code{indSamp=FALSE} in \code{Delta} function call). \cr
 \tab NkMatrix \tab \code{character} \tab Number of samples within each strata (only if \code{indSamp=FALSE} in \code{Delta} function call). \cr
 \tab WkMatrix \tab \code{character} \tab Sampled weight within each strata in grams (only if \code{indSamp=FALSE} in \code{Delta} function call). \cr
 \tab SampDeltaMat \tab \code{character} \tab Data.frame of delta values within each sample (only if \code{indSamp=TRUE} in \code{Delta} function call). \cr
 \tab tab \tab \code{character} \tab Data.frame resulting from treatment (only if \code{indSamp=TRUE} in \code{Delta} function call). \cr
 \tab DFsamp \tab \code{character} \tab Informations about each sample (only if \code{indSamp=TRUE} in \code{Delta} function call).
}
}



\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}
\author{Mathieu Merzereaud}
\seealso{\code{\link{DeltaID}}, \code{\link{plot.Delta}}, \code{\link{plot.DeltaID}}, \code{\link{plot.LD}}, \code{\link{plot.Samp}}
}

\examples{
showClass("DeltA")
}
\keyword{classes}

