\name{plot.ALmult-methods}
\docType{methods}
\alias{plot,ALmult-method}
\alias{plot.ALmult}
\title{Plot "ALmult" object}
\description{
This method plots the output of \emph{AL.multi} procedure as an object of class \emph{ALmult}.
}

\usage{\S4method{plot}{ALmult}(x,show.legend="right",\dots)}

\arguments{
  \item{x}{An \emph{ALmult} object created by \emph{AL.multi} procedure.}
  \item{show.legend}{Display the legend (\code{""} means "no legend").}
  \item{...}{Further graphical parameters.}
}

\section{Methods}{
\describe{
	\item{plot}{\code{signature(x="ANY",y="ANY")}: generic function: see 'plot'.}
	\item{plot}{\code{signature(ALmult)}: plotting procedure of an object of class \emph{ALmult}.}
}}


\references{Gerritsen, H.D., McGrath, D., and Lordan, C. (2006)
\emph{A simple method for comparing age-length keys reveals significant regional 
differences within a single stock of haddock (Melanogrammus aeglefinus)}. ICES Journal of Marine Science, 63: 1096-1100.
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{ALmult}}, \code{\link{AL.multi}}, \code{\link{GraphsPar}}, \code{\link{plot}}
}

\examples{
data(sole3.cs)
obj <- AL.multi(sole3.cs,tempStrata="quarter",spaceStrata="rect",
                elmts=list(tp=c("2","3"),sp="all",tc="all"),
                grps="tempStrata",age.plus=6)
plot(obj,l.col=c("steelblue","violetred2"))
}
\keyword{dplot}
