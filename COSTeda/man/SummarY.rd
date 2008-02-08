\name{SummarY-methods}
\docType{methods}
\alias{SummarY}
\alias{SummarY,csData-method}
\alias{SummarY,clData-method}
\alias{SummarY,ceData-method}
\alias{SummarY-methods}
\title{Summary "plus" procedure for "csData", "clData" and "ceData" objects.}
\description{
These methods implements a special \emph{summary} procedure for objects of class \emph{csData}, \emph{clData} and \emph{ceData}.
}

\usage{
\S4method{SummarY}{csData}(object,tab="tr",sizeMax=20,except=NULL,\dots)
\S4method{SummarY}{clData}(object,sizeMax=20,except=NULL,\dots)
\S4method{SummarY}{ceData}(object,sizeMax=20,except=NULL,\dots)
}

\arguments{
  \item{object}{An object of class \emph{csData}, \emph{clData} or \emph{ceData}.}
  \item{tab}{Character specifying one slot of \emph{csData} object.}
  \item{sizeMax}{Numerical value specifying the number of rows of output.}
  \item{except}{Character specifying fields to omit.}
  \item{...}{Further parameters.}
}

\section{Methods}{
\describe{
	\item{SummarY}{\code{signature(csData)}: summary for a \emph{csData} object.}
	\item{SummarY}{\code{signature(clData)}: summary for a \emph{clData} object.}
  \item{SummarY}{\code{signature(ceData)}: summary for a \emph{ceData} object.}
}}


\author{Mathieu Merzereaud}

\seealso{\code{\link{summary}}
}

\examples{
data(sole3.cs)
SummarY(sole3.cs,tab="ca")
SummarY(sole3.cs,tab="ca",sizeMax=28)
SummarY(sole3.cs,tab="ca",except="trpCode")
}

\keyword{manip}