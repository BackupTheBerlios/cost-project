\name{WL.plot,csData-method}
\alias{WL.plot}
\alias{WL.plot,csData-method}
\docType{methods}
\title{Plots of individual biological parameters : Weight-Length}
\description{
This method implements a scatter plot of individual weight-length data. It requires a \emph{csData} object with \emph{ca} table.
It also allows outliers identification.
}

\usage{
WL.plot(object,selection=FALSE,\dots)
}

\arguments{
  \item{object}{A \emph{csData} object with \emph{ca} informations.}
  \item{selection}{If \code{TRUE}, outlier(s) identification can be done and corresponding data is returned.}
  \item{...}{Further graphical arguments}
}

\author{Mathieu Merzereaud}

\seealso{\code{\link{ML.plot}}, \code{\link{ML.boxplot}}, \code{\link{SL.boxplot}}, \code{\link{L.plot.design}}
}

\examples{
data(sole)
WL.plot(sole.cs)

#4 graphs on the same page
library(grid)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))

pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))
print(WL.plot(sole.cs),newpage=FALSE) ; popViewport(1)
pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
print(ML.plot(sole.cs),newpage=FALSE) ; popViewport(1)
pushViewport(viewport(layout.pos.col=1, layout.pos.row=2))
print(SL.boxplot(sole.cs),newpage=FALSE) ; popViewport(1)
pushViewport(viewport(layout.pos.col=2, layout.pos.row=2))
print(ML.boxplot(sole.cs),newpage=FALSE) ; popViewport(1)
}
\keyword{methods}
