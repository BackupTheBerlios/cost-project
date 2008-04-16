\name{dltCalc,csData-method}
\alias{dltCalc}
\alias{dltCalc,csData-method}
\docType{methods}
\title{Calculation of Delta index for sampling outliers detection and variance calculation}
\description{
This method implements the calculation of Delta values, derived from the formulation of the variance in landings-at-length.
It requires a \emph{csData} object built from \pkg{COSTcore} package. Length distribution informations are taken from \emph{hl} table.
}

\usage{
dltCalc(object,species,timeStrata=NULL,spaceStrata=NULL,techStrata=NULL,
      indSamp=FALSE,strategy="metier",fraction="LAN",\dots)
}

\arguments{
  \item{object}{A \emph{csData} object with \emph{tr}, \emph{hh}, \emph{sl} and \emph{hl} informations.}
  \item{species}{Field specifying species (e.g \code{"Solea vulgaris"}).}
  \item{timeStrata}{Field specifying time stratification (e.g \code{"quarter"}).}
  \item{spaceStrata}{Field specifying space stratification (e.g \code{"area"}).}
  \item{techStrata}{Field specifying technical stratification (e.g \code{"commCat"}).}
  \item{indSamp}{If \code{TRUE}, output is within each sample and is dedicated to outliers detection. If \code{FALSE}, output is within length classes and is dedicated to variance calculation.}
  \item{strategy}{To be chosen between \code{"metier"} and \code{"cc"} (for commercial categories).}
  \item{fraction}{Fate of the catch on which calculation is made. To be chosen between \code{"LAN"}, \code{"DIS"} and \code{"all"}.}
  \item{...}{Further arguments.}
}

\details{For more informations about arguments, see \emph{FishFrame/COST Exchange format specification}.}

\value{An object of class \emph{dltCls}.}

\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}
\author{Mathieu Merzereaud}
\seealso{\code{\link{dltCls}}, \code{\link{dltId}}, \code{\link{dltPlot}}, \code{\link{plot.dltId}}, \code{\link{smpPlot}}, \code{\link{lenDisPlot}}
}

\examples{
data(sole)
obj <- dltCalc(sole.cs,"Solea vulgaris",timeStrata="quarter",techStrata="commCat",
               indSamp=TRUE,strategy="cc")  
}
\keyword{methods}
