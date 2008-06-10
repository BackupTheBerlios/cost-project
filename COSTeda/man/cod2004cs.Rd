\name{cod2004cs}
\alias{cod2004cs}
\alias{had2004cs}
\alias{sai2004cs}
\alias{whi2004cs}
\docType{data}
\title{FRS commercial sampling data }
\description{
FRS commercial sampling data in the COST data exchange format.
\cr Currently consists of cod, haddock, whiting and saithe for 2004. 
}
\usage{
data(cod2004cs)
data(had2004cs)
data(whi2004cs)
data(sai2004cs)
}
\format{
Formal class 'csData' [package "COSTcore"] objects with 6 slots
\cr@desc: description
\cr@tr: data.frame of 15 variables 
\cr@hh: data.frame of 29 variables 
\cr@sl: data.frame of 16 variables
\cr@hl: data.frame of 16 variables
\cr@ca: data.frame of 28 variables
\cr see csData for details of the variables   
}
\details{
Commercial sampling data from vessels landing the names species into Scottish ports during 2004,  
aggregated by month, ICES statistical rectangle, ICES area and gear type.  
Compiled from FRS sampling programme data. 
}
\section{Warning }{These data sets are for test purposes only and should not to be used or cited without prior permission.}

\source{
  FRS Marine Laboratory, Aberdeen. Scotland. 
}



%\references{
%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(cod2004cs)
dim(cod2004cs)
str(cod2004cs)
%## maybe str(cod2004cs) ; plot(cod2004cs) ...
}
\keyword{datasets}
