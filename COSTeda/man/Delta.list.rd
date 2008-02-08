\name{Delta.list-class}
\docType{class}
\alias{Delta.list}
\alias{Delta.list-class}
\title{Class "Delta.list"}
\description{Any object of this class is created by \emph{Delta.cs} procedure.}

\section{Slots}{
\tabular{lrll}{
\bold{slot} \tab  \tab \bold{class} \tab \bold{description} \cr
\bold{\code{Dlist}} \tab \tab \code{list} \tab Stratification, length-weight distributions and Delta \cr
                                \tab \tab \tab values calculated by \emph{Delta.cs} procedure. \cr
 \tab species \tab \code{character} \tab Specified species. \cr
 \tab fraction \tab \code{character} \tab Specified catch category. \cr
 \tab strata \tab \code{list} \tab Specified stratification with components. \cr 
 \tab lenCode \tab \code{character} \tab Length classes step. \cr
 \tab d.jku \tab \code{matrix} \tab Length distributions by FO. \cr
 \tab w.ku \tab \code{numeric} \tab Fraction weight by FO. \cr
 \tab delta.jku \tab \code{matrix} \tab Delta values by length class and FO. \cr
 \tab Delta \tab \code{data.frame} \tab Delta values by FO with corresponding stratas values. 
}}

\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "Delta.list")}: graphical display of Delta values by FO for specified stratification, and outlier identification.}
    \item{plot.Delta.list}{\code{signature(x = "Delta.list")}: See above.}
    \item{plot.Samp}{\code{signature(x = "Delta.list")}: graphical display of length distribution for one given sample, compared with overall relative length distribution (individualized \emph{plot.Delta.length} procedure).}
	 }
}


\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}
\author{Mathieu Merzereaud}
\seealso{\code{\link{Delta.length}}, \code{\link{plot.Delta.list}}, \code{\link{plot.Delta.length}}, \code{\link{Delta.cs}}, \code{\link{plot.LD}}, \code{\link{plot.Samp}}
}

\examples{
showClass("Delta.list")
}
\keyword{classes}

