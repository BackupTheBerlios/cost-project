\name{alkLgthRecSim}
\alias{alkLgthRecSim}
\docType{methods}
\title{Method for managing gaps in age-length keys for simulated data sets}
\description{
This function is the counterpart of the \code{alkLgthRec} function 
for \code{simDataCons} objects. It provides various methods to solve alk gaps problems. 
Input object is updated according to chosen method(s) (grouped/recoded length classes, addition of 'virtual' individuals,...)  
}

\usage{
alkLgthRecSim(object,type="stepIncr",value,preview=FALSE,postview=TRUE,update=FALSE,\dots)
}
                           
\arguments{
  \item{object}{A \code{simDataCons} object with simulated data sets}
  \item{type}{Character for chosen method. Values are :
    \item{"stepIncr"}{Default parameter. Length class step is increased to specified \code{value} parameter (default value=10)}
    \item{"fillMiss"}{All gaps (with size <= value) are filled out with the sum of surrounding recorded classes (default value=1)}
%    \item{"sExtrGrp"}{The 'value' first classes are grouped (default value=1)}
%    \item{"lExtrGrp"}{The 'value' last classes are grouped (default value=1)}
    \item{"sFillMiss"}{The 'value' empty classe(s) prior to first recorded length class is filled out with the latter (default value=1)}
    \item{"lFillMiss"}{The 'value' empty classe(s) following last recorded length class is filled out with the latter (default value=1)} 
%    \item{"sFillAge"}{The 'value' empty classe(s) prior to first recorded length class is filled out with 1 individual of minimal age (default value=1)}
} 
  \item{value}{Numerical parameter for chosen method (see 'type').}
  \item{preview}{Logical. If \code{TRUE}, original age length key is displayed.}
  \item{postview}{Logical. If \code{TRUE}, new age length key is displayed.}
  \item{update}{Logical. If \code{TRUE}, 'csDataCons' object is updated in accordance with chosen method, and then returned. 
  If \code{FALSE}, descriptive elements about updated alk are returned (see 'values'), but input object remains unchanged.}
  \item{...}{Further arguments, and particularly a \code{start} numerical parameter specifying the first considered length class when recoding (only useful for 'type="stepIncr"'). 
  Default value is the minimum aged length class in \code{ca} table.}
}

\value{If \code{update=FALSE}, returned elements within \code{@samples} of the \code{simDataCons} object are : \code{$alk} is the raw resulting age-length key, \code{$propMiss} are short statistics about gaps (see 'propMissLgthCons' method), 
\code{$lgthCls} is a description of length classes recoding for 'stepIncr', 'sExtrGrp' and 'lExtrGrp' methods and \code{$addIndTab} is a description of added virtual individuals 
for other methods.
}

\author{Dorleta garcia \email{dgarcia@azti.es}}        

\seealso{
\code{\link[COSTdbe]{alkLgthRec}}
}

%\examples{ }

\keyword{methods}
