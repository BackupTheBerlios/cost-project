\name{tabConsist}
\alias{tabConsist}

\title{Field composition over COST objects}

\description{
This function gives a description of a specified field contents in various COST objects. 
}

\usage{
tabConsist(lTab,field)
}

\arguments{
  \item{lTab}{ A list of COST objects (within the same class).}
  \item{field}{ Character specifying the field to describe.}
}


\author{Mathieu Merzereaud}

\seealso{\code{\link{dfApply}}}
\examples{

data(sole)
tabConsist(list(sole.cs,sole.ce,sole.cl),"area")

}

\keyword{manip}

