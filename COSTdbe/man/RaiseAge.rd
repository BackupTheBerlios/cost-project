\name{RaiseAge}
\alias{RaiseAge}
\alias{RaiseAge,dbeOutput,csDataCons-method}
\docType{methods}
\title{Estimation of total numbers-at-age from market sampling}
\description{
ToDo
}

\usage{
RaiseAge(dbeOutput,csObject,type="fixed",sex=as.character(NA),\dots)
}

\arguments{
  \item{dbeOutput}{A \emph{dbeOutput} object.}
  \item{csObject}{A \emph{csDataCons} object matching 'dbeOutput' specifications.}
  \item{type}{Allocation strategy used to establish the age-length data.}
  \item{sex}{Sex}
  \item{\dots}{Further arguments}  
}

\references{ToDo}

\value{An updated object of class \emph{dbeOutput}.}


\author{Mathieu Merzereaud}
\seealso{\code{\link{dbeOutput}}, \code{\link{dbeObject}}, \code{\link[COSTcore]{csDataCons}}, \code{\link[COSTcore]{clDataCons}}
}

\examples{

}
\keyword{methods}
