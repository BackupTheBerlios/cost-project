\name{dltPlot,csData-method}
\docType{methods}
\alias{dltPlot}
\alias{dltPlot,csData-method}
\title{Delta values plotting procedure}
\description{
Plotting procedure of delta values calculated with a call to \code{dltCalc} method. It allows outliers identification to create a new object of class \emph{dltId}.
}

\usage{dltPlot(x,species,timeStrata=NULL,spaceStrata=NULL,techStrata=NULL, 
        elmts=list(tp="all",sp="all",tc="all"),strategy="metier",fraction="LAN",
        strat1,strat2="NULL",selection=FALSE,show.legend="right",shift=FALSE,\dots)}

\arguments{
  \item{x}{A \emph{csData} object.}
  \item{species}{Field specifying species (e.g \code{"Solea vulgaris"}).}
  \item{timeStrata}{Field specifying time stratification (e.g \code{"quarter"}).}
  \item{spaceStrata}{Field specifying space stratification (e.g \code{"area"}).}
  \item{techStrata}{Field specifying technical stratification (e.g \code{"commCat"}).}
  \item{elmts}{List of levels for specified stratification to be displayed on the graph.}
  \item{strategy}{To be chosen between \code{"metier"} (default value) and \code{"cc"}.}
  \item{fraction}{Fate of the catch on which calculation is made. To be chosen between \code{"LAN"}, \code{"DIS"} and \code{"all"}.}
  \item{strat1}{Optionnal. To be chosen between \code{"timeStrata"}, \code{"spaceStrata"} and \code{"techStrata"}. Primary stratification for graphical display.}
  \item{strat2}{Optionnal. To be chosen between \code{"timeStrata"}, \code{"spaceStrata"} and \code{"techStrata"}. Secondary stratification for graphical display.}   
  \item{selection}{If \code{TRUE}, outliers identification is made, and a \emph{dltId} object is returned. 
  Displayed values during identification process are taken from 'SampNum' field from \code{sampId} returned data.frame.}
  \item{show.legend}{Display the legend (\code{""} means "no legend").}
  \item{shift}{If \code{TRUE}, displayed text is shifted.}
  \item{...}{Further graphical parameters.}
}



\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}

\author{Mathieu Merzereaud}
\seealso{\code{\link{dltCls}}, \code{\link{dltId}}, \code{\link{dltCalc}}, \code{\link{plot.dltId}}, \code{\link{smpPlot}}, \code{\link{lenDisPlot}}
}

\examples{
data(sole)
dltPlot(sole.cs,species="Solea vulgaris",timeStrata="quarter",techStrata="commCat",
        strat1="techStrata",strat2="timeStrata",strategy="cc")
}

\keyword{dplot}
