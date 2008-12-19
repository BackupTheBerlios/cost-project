\name{CEFAS2006cs}
\alias{CEFAS2006cs}
\alias{CEFAS2006ce}
\alias{CEFAS2006cl}
\docType{data}
\title{CEFAS sampling, landings and effort data}
\description{
CEFAS commercial sampling, landings and effort data in the COST data exchange format.
\cr Mixed spaecies data for 2006. 
}
\usage{
data(CEFAS2006cs)
data(CEFAS2006ce)
data(CEFAS2006cl)
}
\format{

Formal class 'csData' [package "COSTcore"] with 6 slots: \code{desc}, \code{tr}, \code{hh}, \code{sl}, \code{hl} and \code{ca}.

Formal class 'ceData' [package "COSTcore"] with 2 slots: \code{desc} and \code{ce}.

Formal class 'clData' [package "COSTcore"] with 2 slots: \code{desc} and \code{cl}.

}
\details{
The tables below give details on the fields of the Cefas Biological Sampling System (BSS) and Fishing Activity Database (FAD) used in producing market sampling data for the FishFrame-COST data exchange format description 26Feb2008.

Note that:
 data for some of the mandatory fields are not stored so may have to be assumed or  imputed from other fields
species and gear are currently in local codes and landings value is in pounds sterling, we will look to convert these in a later extraction or within R.

David Maxwell and Peter Robinson, 2 April 2008.

}
\section{Warning }{
Cefas data provided for COST (EU FISH/2006/15 Lot 2)
The data are not official statistics. Please contact the Marine Fisheries Agency if official catch or landings statistics are required.
The data are not to be used outside the current project (EU FISH/2006/15 Lot 2) or stored after the project.
Defra must be acknowledged when the data are used.
Cefas must give permission for any publication involving the data. 

}

\source{
  CEFAS 
}
%\references{
%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(CEFAS2006cs)

}
\keyword{datasets}
