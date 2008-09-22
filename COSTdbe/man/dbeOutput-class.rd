\name{dbeOutput-class}
\docType{class}
\alias{dbeOutput}
\alias{dbeOutput-class}
\title{Class "dbeOutput"}
\description{Outcome object from \emph{COSTdbe} methods}

\section{Slots}{
\tabular{lllll}{
\bold{slot} \tab \bold{desc} \tab \bold{elements} \tab \bold{class} \tab \bold{description} \cr
\bold{\code{desc}} \tab \tab \tab \code{character} \tab Descriptive slot \cr
\bold{\code{species}} \tab \tab \tab \code{character} \tab Species \cr
\bold{\code{catchCat}} \tab \tab \tab \code{character} \tab Catch category (eg "LAN", "DIS" or "all") \cr
\bold{\code{param}} \tab \tab \tab \code{character} \tab Parameter estimated (eg "maturity","sex-ratio",...)  \cr
\bold{\code{strataDesc}} \tab \tab \tab \code{strIni} \tab Stratification considered \cr
\bold{\code{methodDesc}} \tab \tab \tab \code{character} \tab Used method (eg "analytical","bootstrap" or "Bayesian") \cr
\bold{\code{desc}} \tab \tab \tab \code{character} \tab Descriptive slot \cr
}
All the following slots contain dataframes with at least 4 fields : \emph{time}, \emph{space}, \emph{technical} and \emph{value}. \cr
\tabular{lllll}{
\bold{nSamp} \tab \tab \tab \code{data.frame} \tab Number of samples\cr
\bold{nMes} \tab \tab \tab \code{data.frame} \tab Number of individual measured\cr
\bold{lenStruc} \tab \tab \tab \code{list} \tab Estimates of the length structure\cr
 \tab \tab estim \tab \code{data.frame} \tab Final estimates (\emph{length} field added). \cr
 \tab \tab rep \tab \code{data.frame} \tab Resampling replicates (\emph{length} and \emph{iter} fields added). \cr
\bold{lenVar} \tab \tab \tab \code{data.frame} \tab Estimates of the variance of 'lenStruc'\cr
\bold{ageStruc} \tab \tab \tab \code{list} \tab Estimates of the age structure\cr
 \tab \tab estim \tab \code{data.frame} \tab Final estimates (\emph{age} field added). \cr
 \tab \tab rep \tab \code{data.frame} \tab Resampling replicates (\emph{age} and \emph{iter} fields added). \cr
\bold{ageVar} \tab \tab \tab \code{data.frame} \tab Estimates of the variance of 'ageStruc'\cr
\bold{totalN} \tab \tab \tab \code{list} \tab Estimates of the total number of the parameters\cr
 \tab \tab estim \tab \code{data.frame} \tab Final estimates. \cr
 \tab \tab rep \tab \code{data.frame} \tab Resampling replicates (\emph{iter} field added). \cr
\bold{totalNvar} \tab \tab \tab \code{data.frame} \tab Estimates of the variance of 'totalN'\cr
\bold{totalW} \tab \tab \tab \code{list} \tab Estimates of the total weight of the parameters\cr
 \tab \tab estim \tab \code{data.frame} \tab Final estimates. \cr
 \tab \tab rep \tab \code{data.frame} \tab Resampling replicates (\emph{iter} field added). \cr
\bold{totalWvar} \tab \tab \tab \code{data.frame} \tab Estimates of the variance of 'totalW'\cr
}}

\author{Mathieu Merzereaud}

\seealso{\code{\link{totVolume}}, \code{\link{dbePlot}}
}

\examples{
showClass("dbeOutput")
}

\keyword{classes}

