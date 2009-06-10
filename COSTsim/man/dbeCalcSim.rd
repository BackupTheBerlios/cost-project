\name{dbeCalcSim}
\alias{dbeCalcSim}
\alias{dbeCalcSim,dbeOutputSim-method}
\docType{methods}                       
\title{CI and CV calculation} 
\description{
Method for calculating coefficients of variation or confidence intervals from dbeOutputSim object estimates. This function is equivalent
to \code{dbeCalc} for dbeOutput but with the \code{update} argument always set equal to TRUE
}

\usage{
dbeCalcSim(object,type="CI",vrbl="l",probs=c(0.025,0.975),replicates=FALSE,\dots)
}

\arguments{
  \item{object}{A \code{dbeOutputSim} object.}
  \item{type}{Character. "CI" for confidence interval calculation, or "CV" for coefficients of variation.}
  \item{vrbl}{Character specifying 'dbeOutput' estimates on which calculation is applied : "l" for length structure, "a" for age structure, "n" for total number estimates, "w" for total weight estimates.}
  \item{probs}{Numeric vector of probabilities with values in [0,1]. Defines CI bounds (relevant only if \code{type="CI"}). See \code{quantile}.}
  \item{replicates}{Logical. If \code{TRUE}, calculation is made from @...\$rep elements ; if \code{FALSE}, @...\$estim and @...Var 'dbeOutput' data are used.}
  \item{...}{Further arguments used as '\code{quantile} method input parameter (if \code{type="CI"} and besides \code{probs} parameter).}
}

\details{
If calculation is made from replicates (see \emph{replicates} parameter), confidence interval is estimated using \emph{quantile} function with \emph{probs} and \dots parameters.
If calculation is made from estimates, normal distribution of total estimates is assumed to compute confidence intervals. 
Possible resulting negative bounds are automatically replaced by 0 in output object. 
}

\author{Dorleta Garcia \email{dgarcia@azti.es}}

\seealso{\code{\link{dbeOutputSim}}, \code{\link[COSTdbe]{dbeCalc}} \code{\link{quantile}}
}

%\examples{}

\keyword{methods}
