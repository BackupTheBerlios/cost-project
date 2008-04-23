\name{space.plot}
\alias{space.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial plotting function}
\description{
This function maps a numeric variable grouped by a function over spatial strata. The possible spatial strata are ICES 
statistical rectangle or ICES (FAO) areas, covering a geographic area extending west to -50 degrees longitude and 
east to 70 degrees, and between 20 and 85 degrees of latitude. This function is based on \code{\link{tapply}} 
and is designed for plotting COST objects but is capable of stand alone use.   
}
\usage{
space.plot(variable, SpaceStrat, func, xlim = NULL, ylim = NULL, zlim = NULL, 
xlab = NULL, ylab = NULL, breaks = NULL, maptype = "image", plotmap = TRUE, 
overlay = FALSE, ices.divs = FALSE, depths = FALSE, statrects = FALSE, 
fcoast = FALSE, landmass = FALSE, pch = 1, colour = TRUE, col.coast = "blue", 
col.cont = "grey", col.pch = "red", col.rect = "grey", col.land = "snow2", 
col.depth = "grey", col.text = 1, scale = FALSE, scale.title = "", 
scale.cex = 0.6, scaleplace = "bottomright", scale.box="o",cex.max.bubble = 2, 
threshold = 0, digits.text = 0, cex.text = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variable}{ A vector of values of the numeric variable to be plotted.}
  \item{SpaceStrat}{A vector of the spatial strata over which \code{variable} is to be grouped, of the same length as \code{variable}. 
  Options are ICES statistical rectangle in '45D3' format, or ICES area in Roman numeric format 'VIIa'.}
  \item{func}{ The grouping function to be applied to \code{variable} over \code{SpaceStrat}. This can be any 
  function that can be passed as the \code{func} argument to tapply.
  }
  \item{xlim}{Vector of c(min, max) for the longitude limits of the plotting area in decimal degrees. 
  Longitudes to the west of 0 degrees are expressed as negatives.}
  \item{ylim}{Vector of c(min, max) for the latitude limits of the plotting area in decimal degrees. }
  \item{zlim}{Vector of c(min, max) for the limits to \code{variable}}
  \item{xlab}{Label for the x axis}
  \item{ylab}{Label for the y axis}
  \item{breaks}{The breaks by which to split \code{variable}. The default is NULL which results in seq(zmin, zmax, length.out=8) being used.}
  \item{maptype}{One of "image", "contour", "bubble" or "values". The default is "image". 
  \cr "image" plots \code{variable} grouped by \code{func} as an image plot over the chosen \code{SpaceStrat}.
  \cr "contour" plots \code{variable} grouped by \code{func} as contour lines over the chosen \code{SpaceStrat}. 
  \cr "bubble" plots \code{variable} grouped by \code{func} as a scaled plotting character at the (approximate) 
  centre of the chosen \code{SpaceStrat}.
  \cr "values" plots \code{variable} grouped by \code{func} as numeric values at the (approximate) centre of the chosen \code{SpaceStrat}.}
  
  \item{plotmap}{Logical. If FALSE plots the map area but omits the "image", "contour", "bubble", or "values".}
  \item{overlay}{Logical. If TRUE adds to the proceeding plot.}
  \item{ices.divs}{Logical. If TRUE adds ICES division lines.}
  \item{depths}{Vector of the form c(50,100,200,500). Adds the specified depth contour.}
  \item{statrects}{Logical. If TRUE adds ICES statistical rectangle division lines.}
  \item{fcoast}{Logical. If TRUE adds 1:250,000 scale coastline. Not available when landmass=TRUE.}
  \item{landmass}{Logical. If TRUE adds landmass polygons.}
  \item{pch}{Plotting character used for maptype="bubble", the default is pch=1.}
  \item{colour}{Logical. If TRUE the colours used for image plots are rev(heatcolors(n-1)) where n is the 
  length of breaks. If FALSE uses grey(rev(seq(0.1, 0.9, length.out =(n-1)))). }
  \item{col.coast}{Colour to be used for plotting the coastline, the default is "blue".}
  \item{col.cont}{Colour to be used for the contour lines, the default is "grey".}
  \item{col.pch}{Colour to be used for the plotting character, the default is "red".}
  \item{col.rect}{Colour to be used for statistical rectangle lines, the default is "grey".}
  \item{col.land}{Colour to be used for the landmass polygons, the default is "snow2".}
  \item{col.depth}{Colour to be used for the depth contour lines, the default is "grey".}
  \item{col.text}{Colour to be used for values when maptype="values", the default is 1 (black).}
  \item{scale}{Logical. If TRUE a scale is added to the plot.}
  \item{scale.title}{Title for the scale, the default is no title.}
  \item{scale.cex}{Cex to be used for the size of the scale box.}
  \item{scaleplace}{Location of the scale, default is "bottomright", available options are "bottomright", "bottom", "bottomleft", "left",
   "topleft", "top", "topright", "right" and "center".}
  \item{scale.box}{The type of box to be drawn around the scale, "o" the default or "n" none.}
  \item{cex.max.bubble}{Maximum size of the plotting character when maptype="bubble", the default is 2.}
  \item{threshold}{Only relevant when maptype="values". The value below which points are shown with a "+" 
  symbol rather than their actual value. The default is zero.}
  \item{digits.text}{Number of decimal places to print when maptype="values", the default is 0.}
  \item{cex.text}{Cex to be used for printed values when maptype="values". The default is 1.}
  \item{\dots}{Additional arguments such as those to \code{plot}}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Invisibly returns various parameters used or defined as a list: 
  \item{lon}{Longitude values used in the plot, i.e. the x range of the image plot}
  \item{lat}{Latitude values used in the plot, i.e. the y range of the image plot}
  \item{values}{The plotted values after the application of \code{func}, i.e. the z range of the image plot.}
  \item{bks}{The breaks used.}
  \item{cols}{The colours used.}
  \item{maptype}{The maptype used.}
  \item{scale.cex}{The scale.cex used.}  
  \item{cex.max.bubble}{The cex.max.bubble used.}  

}
%\references{ ~put references to the literature/web site here ~ }
\author{Alastair Pout  \email{a.pout@marlab.ac.uk}, Liz Clarke  \email{e.d.clarke@marlab.ac.uk}}
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{strata.space.plot}} for the plotting of COST objects by spatial strata}
\examples{

data(cod2004ce)
# Effort by Scottish vessels for cod in 2004: 
#   number of trips by ICES statistical rectangle  
space.plot(cod2004ce@ce$trpNum,cod2004ce@ce$rect,sum,statrects=TRUE,
scale.title="Number of trips")
title("Effort by Scottish vessels for cod in 2004")
# Landings by Scottish vessels of cod in 2004: 
#   landed weight by ICES statistical rectangle  
data(cod2004cl)
space.plot(cod2004cl@cl$landWt,cod2004cl@cl$rect,sum,statrects=TRUE,
scale.title="Landed weight (kg)")
title("Landings by Scottish vessels of cod in 2004") 



# A comparason of landed weight and number of fish measured by ICES rectangle
data(cod2004cl)
data(cod2004cs)
# adding the date and area fields to the csData tables 
newcs <-mergecsData(cod2004cs)
# first plotting out landings to get the full spatial extent.
space.plot(cod2004cl@cl$landWt,cod2004cl@cl$rect,sum,
maptype="bubble",col.pch=1,cex.max.bubble=4
,scaleplace="topleft",scale=TRUE,scale.box="n",scale.title="Landed wgt (kg)")

# then overlay numbers of measured fish
space.plot(newcs@hl$lenNum,newcs@hl$rect,sum,maptype="bubble",pch=16,
overlay=TRUE,cex.max.bubble=4,scaleplace="bottomleft",
scale=TRUE,scale.box="n",scale.title="Measured fish")

# redoing the landings plot to show landings < sampling  
space.plot(cod2004cl@cl$landWt,cod2004cl@cl$rect,sum,maptype="bubble",
col.pch=1,cex.max.bubble=4,scaleplace="topleft",scale=TRUE,
overlay=TRUE,scale.box="n")

title(sub="Landed weight and number of fish measured, cod 2004") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
