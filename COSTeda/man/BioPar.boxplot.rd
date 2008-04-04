\name{BioPar.boxplot,csData-method}
\alias{BioPar.boxplot}
\alias{BioPar.boxplot,csData-method}
\docType{methods}
\title{Boxplot of individual biological parameters}
\description{
This method implements a boxplot of individual biological data. It requires a \emph{csData} object with \emph{ca} table.
}

\usage{
BioPar.boxplot(object,type="WL",\dots)
}

\arguments{
  \item{object}{A \emph{csData} object with \emph{ca} informations.}
  \item{type}{To be chosen between \code{"WL"} (weight~length), \code{"ML"} (maturity~length) or \code{"SL"} (sex~length).}
  \item{...}{Further graphical arguments.}
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{BioPar.plot}}, \code{\link{L.plot.design}}
}

\examples{
data(sole)
BioPar.boxplot(sole.cs)

#4 graphs on the same page
library(grid)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))

pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))
print(BioPar.plot(sole.cs),newpage=FALSE) ; popViewport(1)
pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
print(BioPar.plot(sole.cs,type="ML"),newpage=FALSE) ; popViewport(1)
pushViewport(viewport(layout.pos.col=1, layout.pos.row=2))
print(BioPar.boxplot(sole.cs,type="SL"),newpage=FALSE) ; popViewport(1)
pushViewport(viewport(layout.pos.col=2, layout.pos.row=2))
print(BioPar.boxplot(sole.cs,type="ML"),newpage=FALSE) ; popViewport(1)

}
\keyword{methods}
