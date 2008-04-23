\name{strata.space.plot}
\alias{strata.space.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Stratified spatial plotting function}
\description{
 A stratified spatial plotting function to plot COST cl ce and cs objects. 
  \cr The possible spatial strata are ICES statistical rectangle and ICES (FAO) areas. 
  \cr The temporal strata are year, quarter and month. 
  \cr The technical strata are national gear codes or level 5 or 6 of the Nantes matrix. 
 }
\usage{
strata.space.plot(object, variable, SpaceStrat, TimeStrat, TechStrat, 
func, nplots=1,multiscale.title,multiscale.cex=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ COST object either a \code{csData}, \code{ceData} or \code{clData} object}
  \item{variable}{ Any numeric variable in the object specified }
  \item{SpaceStrat}{ The spatial strata over which the variable is to be plotted, one of either ICES statistical 
  rectangle "rect", or ICES (FAO) area "area".}
  \item{TimeStrat}{ The temporal stratification to be used, options are all, which is the default, year "year", quarter "quarter" or month "month".}
  \item{TechStrat}{ The technical stratification to be used, options are all, which is the default, "foCatNat", "foCatEu5", and "foCatEu6".  }
  \item{func}{ The grouping function to be applied to \code{variable} over the \code{SpaceStrat}, \code{TimeStrat} \code{TechStrat} 
  combination specified. This can be any function that can be passed as the \code{func} argument to \code{\link{tapply}}}
  \item{nplots}{The number of plots per page. Available multiples are 1, 2, 4, 6, 9 or 12. Some multiple stratifications can 
  result in a large number of cells. Specifying the number of plots per page simply allows control over the presentation of these 
  stratifications. The default is \code{nplots}=1. Where \code{TimeStrat} = "quarter" and \code{TechStrat} is unspecified \code{nplots}=4 is 
  the logical choice. Where \code{TimeStrat} = "month" and \code{TechStrat} is unspecified \code{nplots}=12 is the logical choice. This 
  argument determines the preset \code{\link{layout}} setting to be used, an alternative is to specify the desired layout using 
  par(mfrow=c(...,..)) prior to calling \code{strata.space.plot}. }
  \item{multiscale.title}{Title to be added to the scale produced when \code{nplot} is greater than one.}
   \item{multiscale.cex}{Cex size for the scale, the default is 1.}
  \item{\dots}{ Additional arguments such as those to \code{\link{space.plot}}.}
}
\details{
The spatial extent of the plot area is west to -50 degrees longitude and east to 70 degrees. Latitude 20 to 85 degrees. 
}
%\value{
%  ~Describe the value returned
%  If it is a LIST, use
%  \item{comp1 }{Description of 'comp1'}
%  \item{comp2 }{Description of 'comp2'}
%  ...
%}
%\references{ ~put references to the literature/web site here ~ }
\author{Alastair Pout  \email{a.pout@marlab.ac.uk}}
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{space.plot}} for non-stratified plotting of ceData, clData and csData objects}
\examples{
data(cod2004cl)
# A plot of total landed weight by ICES statistical rectangle for Scottish 
# cod in 2004. This replicates the example for space.plot
strata.space.plot(cod2004cl,"landWt","rect",TimeStrat=NULL,TechStrat=NULL,
sum,scale.title="Landed Weight (kg)",statrects=TRUE,scale=TRUE)

# The same plot stratified by quarter
strata.space.plot(cod2004cl,"landWt","rect",TimeStrat="quarter",
TechStrat=NULL,sum,statrects=TRUE,nplots=4,
multiscale.title="Landed weight (kg) of cod by Scottish vessels 2004")


# The number of sampling events is the length of any of the numeric 
# variables within the tr table. 
# A greyscale plot of number of sampling events by area and 
# level 5 Nantes matrix gear type
data(cod2004cs)
strata.space.plot(cod2004cs,"vslId","area",TimeStrat=NULL,
TechStrat="foCatEu5",length,scale=FALSE,nplots=4,colour=FALSE)
title("Number of samples for Scottish cod 2004",line=2)

}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
