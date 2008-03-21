\name{plot.Delta,csData-method}
\docType{methods}
\alias{plot.Delta}
\alias{plot.Delta,csData-method}
\title{Delta values plotting procedure}
\description{
Plotting procedure of delta values calculated with a call to \code{Delta} method. It allows outliers identification to create a new object of class \emph{DeltaID}.
}

\usage{plot.Delta(x,y=NULL,species,tempStrata=NULL,spaceStrata=NULL,techStrata=NULL, 
        elmts=list(tp="all",sp="all",tc="all"),strategy="metier",strat1,strat2="NULL",selection=FALSE,show.legend="right",\dots)}

\arguments{
  \item{x}{A \emph{csData} object.}
  \item{y}{...}
  \item{species}{Field specifying species (e.g \code{"Solea vulgaris"}).}
  \item{tempStrata}{Field specifying time stratification (e.g \code{"quarter"}).}
  \item{spaceStrata}{Field specifying space stratification (e.g \code{"area"}).}
  \item{techStrata}{Field specifying technical stratification (e.g \code{"commCat"}).}
  \item{elmts}{List of components of stratification. Default value is \code{"all"}.}
  \item{strategy}{To be chosen between \code{"metier"} and \code{"cc"}.}
  \item{strat1}{To be chosen between \code{"tempStrata"}, \code{"spaceStrata"} and \code{"techStrata"}. Primary stratification for graphical display.}
  \item{strat2}{To be chosen between \code{"tempStrata"}, \code{"spaceStrata"} and \code{"techStrata"}. Secondary stratification for graphical display.}   
  \item{selection}{If \code{TRUE}, outliers identification is made, and a \emph{DeltaID} object is returned. 
  Displayed values during identification process are Sample field from \code{tab} object.}
  \item{show.legend}{Display the legend (\code{""} means "no legend").}
  \item{...}{Further graphical parameters.}
}



\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}

\author{Mathieu Merzereaud}
\seealso{\code{\link{DeltaID}}, \code{\link{DeltA-class}}, \code{\link{Delta}}, \code{\link{plot.DeltaID}}, \code{\link{plot.LD}}, \code{\link{plot.Samp}}, \code{\link{plot}}
}

\examples{
data(sole)
plot.Delta(sole.cs,species="Solea vulgaris",tempStrata="quarter",techStrata="commCat",strat1="techStrata",strat2="tempStrata",selection=FALSE,show.legend="right")
}

\keyword{dplot}
