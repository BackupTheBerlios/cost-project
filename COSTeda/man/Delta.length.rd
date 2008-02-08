\name{Delta.length-class}
\docType{class}
\alias{Delta.length}
\alias{Delta.length-class}
\title{Class "Delta.length"}
\description{Any object of this class is created by \emph{plot.Delta.length} procedure.}

\section{Slots}{
\tabular{lrll}{
\bold{slot} \tab  \tab \bold{class} \tab \bold{description} \cr
\bold{\code{Dlength}} \tab \tab \code{list} \tab Stratification, length-weight distributions and Delta \cr
                                  \tab \tab \tab values from parent \emph{Delta.list} object, and informations \cr  
                                  \tab \tab \tab about identified outliers. \cr
 \tab species \tab \code{character} \tab Specified species (from \emph{Delta.list} object). \cr
 \tab fraction \tab \code{character} \tab Specified catch category (from \emph{Delta.list} object). \cr
 \tab strata \tab \code{list} \tab Specified stratification with components (from  \cr
                    \tab \tab \tab \emph{Delta.list} object). \cr
 \tab lenCode \tab \code{character} \tab Length classes step. \cr
 \tab d.jku \tab \code{matrix} \tab Length distributions by FO (from \emph{Delta.list} object). \cr
 \tab w.ku \tab \code{numeric} \tab Fraction weight by FO (from \emph{Delta.list} object). \cr
 \tab delta.jku \tab \code{matrix} \tab Delta values by length class and FO (from \emph{Delta.list}  \cr
                         \tab \tab \tab object). \cr
 \tab Delta \tab \code{data.frame} \tab Delta values by FO with corresponding stratas values \cr  
                         \tab \tab \tab (from \emph{Delta.list} object). \cr
 \tab delta.Id \tab \code{matrix} \tab Delta values by length class and by identified \cr
                        \tab \tab \tab FO (outliers). \cr  
 \tab color.code \tab \code{character} \tab Color codes used for second stratification in outliers \cr
                             \tab \tab \tab identification plot. 
}
}

\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "Delta.length")}: graphical display of length distributions for each outlying sample, compared with overall relative length distribution.}
    \item{plot.Delta.length}{\code{signature(x = "Delta.length")}: See above.}
   }
}


\references{Vigneau, J. and Mahevas, S. (2007)
\emph{Detecting sampling outliers and sampling heterogeneity when catch-at-length is estimated using the ratio estimator}. Oxford Journals.
}
\author{Mathieu Merzereaud}
\seealso{\code{\link{Delta.list}}, \code{\link{plot.Delta.list}}, \code{\link{plot.Delta.length}}, \code{\link{Delta.cs}}, \code{\link{plot.LD}}, \code{\link{plot.Samp}}
}

\examples{
showClass("Delta.length")
}
\keyword{classes}

