\name{plot.ageLenMult}
\docType{methods}
\alias{plot,ageLenMult-method}
\alias{plot.ageLenMult,ageLenMult-method}
\alias{plot.ageLenMult}
\title{Plot "ageLenMult" object}
\description{
This method plots the output of \emph{ageLenMultinom} procedure as an object of class \emph{ageLenMult}.
}

\usage{
\S4method{plot}{ageLenMult}(x,y=NULL,show.legend="right",\dots)}

\arguments{
  \item{x}{An \emph{ageLenMult} object created by \emph{ageLenMultinom} procedure.}
  \item{y}{NULL}
  \item{show.legend}{Display the legend (\code{""} means "no legend").}
  \item{...}{Further graphical parameters.}
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

\seealso{\code{\link{ageLenMult}}, \code{\link{ageLenMultinom}}, \code{\link{GraphsPar}}, \code{\link{plot}}
}

\examples{
data(sole)
sole.cs.val <- csDataVal(sole.cs) 
obj <- ageLenMultinom(sole.cs.val,timeStrata="quarter",spaceStrata="rect",
                elmts=list(tp=c("2","3"),sp="all",tc="all"),
                grps="timeStrata",age.plus=6)
plot(obj,l.col=c("steelblue","violetred2"))
}
\keyword{dplot}
