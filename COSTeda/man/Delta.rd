\name{Delta,csData-method}
\alias{Delta}
\alias{Delta,csData-method}
\docType{methods}
\title{Calculation of Delta index for sampling outliers detection and variance calculation}
\description{
This method implements the calculation of Delta values, derived from the formulation of the variance in landings-at-length.
It requires a \emph{csData} object built from \pkg{COSTcore} package. Length distribution informations are taken from \emph{hl} table.
}

\usage{
Delta(object,species,tempStrata=NULL,spaceStrata=NULL,techStrata=NULL,
      indSamp=FALSE,strategy="metier",\dots)
}

\arguments{
  \item{object}{A \emph{csData} object with \emph{tr}, \emph{hh}, \emph{sl} and \emph{hl} informations.}
  \item{species}{Field specifying species (e.g \code{"Solea vulgaris"}). See Details.}
  \item{tempStrata}{Field specifying time stratification (e.g \code{"quarter"}). See Details.}
  \item{spaceStrata}{Field specifying space stratification (e.g \code{"area"}). See Details.}
  \item{techStrata}{Field specifying technical stratification (e.g \code{"commCat"}). See Details.}
  \item{indSamp}{If \code{TRUE}, output is within each sample. If \code{FALSE}, output is within length classes. Both are complementary.}
  \item{strategy}{To be chosen between \code{"metier"} and \code{"cc"}.}
  \item{...}{Further arguments.}
}

\details{For more informations about arguments, see \emph{FishFrame/COST Exchange format specification}.}

\value{An object of class \emph{DeltA}.}

\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}
\author{Mathieu Merzereaud}
\seealso{\code{\link{DeltaID}}, \code{\link{DeltA-class}}, \code{\link{plot.Delta}}, \code{\link{plot.DeltaID}}, \code{\link{plot.LD}}, \code{\link{plot.Samp}}
}

\examples{
data(sole)
obj <- Delta(sole.cs,"Solea vulgaris",tempStrata="quarter",techStrata="commCat",indSamp=TRUE,strategy="cc")  
}
\keyword{methods}
