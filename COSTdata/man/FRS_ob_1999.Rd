\name{FRS_ob_1999}
%\alias{FRS_ob_trips}
\alias{FRS_ob_1999}
\docType{data}
\title{FRS observer data}
\description{
FRS observer data in the COST data exchange format.
\cr Consists of 53 demersal sampling trips from 1999. Discards length distributions are sampled by haul, the landed length distribution is sampled by trip, and the age given length distribution of a sub-sample of the discarded fraction is pooled by trip.  
}
\usage{
data(FRS_ob_1999)
}
\format{
Formal class 'csData' [package "COSTcore"] objects with 6 slots
\cr@desc: description
\cr@tr: data.frame of 16 variables 
\cr@hh: data.frame of 29 variables 
\cr@sl: data.frame of 17 variables
\cr@hl: data.frame of 16 variables
\cr@ca: data.frame of 31 variables
\cr see csData for details of the variables   
}
\details{
The FRS observer sampling protocol is as follows: \cr
From each haul during the trip the discarded fraction of the catch (consisting of unsorted fish of various species) is 
sub-sampled by the observer; two representative baskets of the discarded fish being obtained. The ratio of the discarded weight to sub-sampled weight is 
estimated by the observer and expressed in terms of "baskests". 
The length frequencies are taken from all cod, haddock, whiting and saithe in the sub-sample. 
The otoliths are collected from cod, haddock, whiting and saithe of each length class in the sub-sampled fraction until (usually) 5 otoliths are obtained for each length class for each species.  
Some length classes are not represented and some have fewer than 5 individuals. 
At the end of the trip the landed fraction of the catch, which will have been sorted into commercial size classes, is sub-sampled and the length frequencies recorded for cod, haddock, whiting and saithe. No otoliths are taken of the landed fraction. 

This observer data therefore consists of three components: 
\cr1. The length distribution for a sub-sample of the discarded fraction of the catch by haul.
This is contained in the hl table where \kbd{\$catchCat} = "DIS", trips are identifiable 
by trip code \kbd{\$trpCode}, and the individual hauls by station number \kbd{\$staNum}
\cr2. The length distribution of a sub-sample of the landed fraction pooled by trip.  
This is contained in the hl table where \kbd{\$catchCat} = "LAN", and the station number is \kbd{\$staNum = 999}. 
\cr3. The age given length distribution of a sub-sampled fraction of the discarded catch pooled by trip. 
This is contained in the ca table where \kbd{\$catchCat} = "DIS", and the station number is \kbd{\$staNum} = 999. 

Note that to obtain a catch weight \kbd{\$wt} and sub-sampled weight \kbd{\$subSampWt} by species in the sl table these values have been obtained retrospectively from the species' length frequency distribution using a standard weight length relationship, and the raising factor of number of discarded "baskets" to number sub-sampled baskets. 

}
\section{Warning }{This is a test data set only and should not to be used or cited without prior permission.}
\source{
  FRS Marine Laboratory, Aberdeen, Scotland. 
}

\examples{
data(FRS_ob_1999)

}
\keyword{datasets}
