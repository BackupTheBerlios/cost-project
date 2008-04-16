\name{ageLenMult-class}
\docType{class}
\alias{ageLenMult}
\alias{ageLenMult-class}
\title{Class "ageLenMult"}
\description{This class contains all informations produced by \emph{ageLenMultinom} method}

\section{Slots}{
\tabular{lrll}{
\bold{slot} \tab  \tab \bold{class} \tab \bold{description} \cr
\bold{\code{objMult}} \tab \tab \code{list} \tab informations from multinomial modelisation of \cr 
                                    \tab \tab \tab age-at-length data \cr
 \tab timeStrata \tab \code{character} \tab time stratification field \cr
 \tab spaceStrata \tab \code{character} \tab space stratification field \cr
 \tab techStrata \tab \code{character} \tab technical stratification field \cr
 \tab grps \tab \code{character} \tab grouping field \cr
 \tab Mm \tab \code{multinomial} \tab a \pkg{nnet} package object resulting from \emph{multinom} procedure \cr
 \tab dat \tab \code{data.frame} \tab Age-at-Length data for specified stratification 
}
}

\section{Methods}{
\describe{
	\item{plot}{\code{signature("ageLenMult")}: plotting procedure of an object of class \emph{ageLenMult}.}
}}



\references{Gerritsen, H.D., McGrath, D., and Lordan, C. (2006)
\emph{A simple method for comparing age-length keys reveals significant regional 
differences within a single stock of haddock (Melanogrammus aeglefinus)}. ICES Journal of Marine Science, 63: 1096-1100.
}

\author{Mathieu Merzereaud}
\seealso{\code{\link{ageLenMultinom}}, \code{\link{plot.ageLenMult}}}

\examples{
showClass("ageLenMult")
}
\keyword{classes}

