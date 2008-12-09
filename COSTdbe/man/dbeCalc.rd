\name{dbeCalc}
\alias{dbeCalc}
\alias{dbeCalc,dbeOutput-method}
\docType{methods}                       
\title{CI and CV calculation} 
\description{
Method for calculating coefficients of variation or confidence intervals from 'dbeOutput' object estimates. Input object can optionally be updated with resulting table.
}

\usage{
dbeCalc(object,type="CI",vrbl="l",probs=c(0.025,0.975),replicates=FALSE,update=FALSE,\dots)
}

\arguments{
  \item{object}{A \emph{dbeOutput} object.}
  \item{type}{Character. "CI" for confidence interval calculation, or "CV" for coefficients of variation.}
  \item{vrbl}{Character specifying 'dbeOutput' estimates on which calculation is applied : "l" for length structure, "a" for age structure, "n" for total number estimates, "w" for total weight estimates.}
  \item{probs}{Numeric vector of probabilities with values in [0,1]. Defines CI bounds (relevant only if \code{type="CI"}). See \emph{quantile}.}
  \item{replicates}{Logical. If \code{TRUE}, calculation is made from @...$rep elements ; if \code{FALSE}, @...$estim and @...Var 'dbeOutput' data are used.}
  \item{update}{Logical. If \code{TRUE}, updated 'dbeOutput' object is returned ; if \code{FALSE}, only resulting dataframe is returned}
  \item{...}{Further arguments used as '\emph{quantile} method input parameter (if \code{type="CI"} and besides \code{probs} parameter).}
}

\details{
If calculation is made from replicates (see \emph{replicates} parameter), confidence interval is estimated using \emph{quantile} fonction with \emph{probs} and \dots parameters.
If calculation is made from estimates, normal distribution of total estimates is assumed.
}


\author{Mathieu Merzereaud}
\seealso{\code{\link{dbeOutput}}, \code{\link{dbePlot}}, \code{\link{quantile}}
}

\examples{
}
\keyword{methods}