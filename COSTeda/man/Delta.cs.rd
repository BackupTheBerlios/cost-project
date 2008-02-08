\name{Delta.cs,csData-method}
\alias{Delta.cs}
\alias{Delta.cs,csData-method}
\docType{methods}
\title{Calculation of Delta index for sampling outliers detection}
\description{
This method implements the calculation of Delta values, derived from the formulation of the variance in landings-at-length.
It requires a \emph{csData} object built from \pkg{COSTcore} package. Length distribution informations are taken from \emph{hl} table.
}

\usage{
Delta.cs(object,species,fraction,strata1,strata2=NULL,elmts1="all",
         elmts2="all",show.data=FALSE,\dots)
}

\arguments{
  \item{object}{A \emph{csData} object with \emph{tr}, \emph{hh}, \emph{sl} and \emph{hl} informations.}
  \item{species}{Field specifying species (e.g \code{"SOL"}). See Details.}
  \item{fraction}{Field specifying catch category (e.g \code{"LAN"}). See Details.}
  \item{strata1}{
  Field specifying primary stratification (to be chosen between \code{"year"},
  \code{"quarter"},\code{"month"},\code{"area"},\code{"rect"},\code{"foCatEu5"},\code{"gear"},...). See Details.}
  \item{strata2}{Field specifying secondary stratification. It can be \code{NULL}. See Details.}
  \item{elmts1}{Components of primary stratification. It can be \code{"all"}. See Details.}
  \item{elmts2}{Components of secondary stratification. It can be \code{"all"}. See Details.}
  \item{show.data}{Display (\code{TRUE}) or only assign result (\code{FALSE}).}
  \item{...}{Further arguments.}
}

\details{For more informations about arguments, see \emph{FishFrame/COST Exchange format specification}.}

\value{An object of class \emph{Delta.list}.}

\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}
\author{Mathieu Merzereaud}
\seealso{\code{\link{Delta.list}}, \code{\link{plot.Delta.list}}, \code{\link{Delta.length}}, \code{\link{plot.Delta.length}}, \code{\link{plot.LD}}, \code{\link{plot.Samp}}
}

\examples{
data(sole3.cs)
object <- sole3.cs
#only sea sampling data is kept
object@tr <- object@tr[object@tr$sampType=="S",]
object@hh <- object@hh[object@hh$sampType=="S",]
object@sl <- object@sl[object@sl$sampType=="S",]
object@hl <- object@hl[object@hl$sampType=="S",]

res <- Delta.cs(object,"SOL","LAN","quarter","month")
}
\keyword{methods}
