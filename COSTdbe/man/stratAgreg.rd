\name{stratAggreg}
\alias{stratAggreg}
\alias{stratAggreg,dbeOutput-method}
\docType{methods}
\title{Aggregation of 'dbeOutput' table}
\description{
This method is aggregating (summing) input 'dbeOutput' object tables over chosen stratifications. 
}

\usage{
stratAggreg(object,timeStrata=TRUE,spaceStrata=TRUE,techStrata=FALSE,\dots)
}

\arguments{
  \item{object}{A \emph{dbeOutput} object.}
  \item{timeStrata}{Logical. If \code{TRUE}, aggregation is made over time strata.}
  \item{spaceStrata}{Logical. If \code{TRUE}, aggregation is made over space strata.}
  \item{techStrata}{Logical. If \code{TRUE}, aggregation is made over technical strata.}
  \item{\dots}{Further arguments}  
}

\references{ToDo}

\value{An updated object of class \emph{dbeOutput}.}
\details{
Aggregating process is only applied to these tables : \code{nSamp}, \code{nMeas}, \code{lenVar}, \code{ageVar}, \code{totalNvar}, \code{totalWvar}, and all \code{estim} list elements.
Other tables in output object are empty (so, for example, 'dbeCalc' method must be called to calculate and insert \emph{cv} table in this new object). 
WARNING : summing variances requires some strong probability assumptions within strata, such as uncorrelated variables. In the 'variance at age' case, this can be problematic 
since the same age-length key is used for every technical strata.   
}

\author{Mathieu Merzereaud}
\seealso{\code{\link{dbeOutput}}, \code{\link{dbeObject}}
}

\examples{

}
\keyword{methods}
