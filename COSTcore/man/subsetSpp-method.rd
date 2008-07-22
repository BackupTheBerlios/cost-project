\name{subsetSpp}
\alias{subsetSpp}
\alias{subsetSpp,csData-method}
\alias{subsetSpp,csDataVal-method}
\docType{methods}
\title{Specific subsetting function applying on SL table from COST objects}
\description{
This method implements subsetting for the raw and validated CS objects provided by COSTcore, proceeding specifically on SL table. This subset only impacts on HL table, 
and preserves the other tables.  
}

\usage{
subsetSpp(x,subset,\dots)
}

\arguments{
  \item{x}{A \emph{csData} or \emph{csDataVal} object.}
  \item{subset}{Logical expression specifying the subsetting to be applied on sl table.}
  \item{...}{Further arguments.}
}

\seealso{\code{\link{subset,csData-method}}
}

\keyword{methods}
