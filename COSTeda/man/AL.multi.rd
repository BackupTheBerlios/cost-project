\name{AL.multi,csData-method}
\alias{AL.multi}
\alias{AL.multi,csData-method}
\docType{methods}
\title{Multinomial modelisation applied to fisheries age-at-length data}
\description{
This method implements a multinomial analysis of age-at-length data, for specified time, space and technical stratification. It requires a \emph{csData} object built from \pkg{COSTcore} package. Information is taken from \emph{ca} table.
}

\usage{
AL.multi(object,tempStrata=NULL,spaceStrata=NULL,techStrata=NULL,
         elmts=list(tp="all",sp="all",tc="all"),grps=NULL,age.plus=-1,\dots)
}

\arguments{
  \item{object}{A \emph{csData} object with \emph{ca} informations.}
  \item{tempStrata}{Field specifying time stratification (e.g \code{"year"},\code{"quarter"},\code{"month"} or \code{NULL}).}
  \item{spaceStrata}{Field specifying space stratification (e.g \code{"area"},\code{"rect"} or \code{NULL}).}
  \item{techStrata}{Field specifying technical stratification (e.g \code{"commCat"} or \code{NULL}).}
  \item{elmts}{List of occurrence(s) for each strata. \emph{tp}, \emph{sp} and \emph{tc} are characters. \code{"all"} includes all occurrences for the strata.}
  \item{grps}{Strata to be differentiated in each panel generated by \code{plot.AL.multi} procedure. \code{NULL} means one graph per strata.}
  \item{age.plus}{Threshold for grouping age (numerical value). (-1)-value means no grouping.}
  \item{...}{Further arguments.}
}

\references{Gerritsen, H.D., McGrath, D., and Lordan, C. (2006)
\emph{A simple method for comparing age-length keys reveals significant regional 
differences within a single stock of haddock (Melanogrammus aeglefinus)}. ICES Journal of Marine Science, 63: 1096-1100.
}

\author{Mathieu Merzereaud}
\seealso{\code{\link{ALmult}}, \code{\link{plot.ALmult}}}

\examples{
data(sole)
AL.multi(sole.cs,age.plus=6)
}
\keyword{methods}
