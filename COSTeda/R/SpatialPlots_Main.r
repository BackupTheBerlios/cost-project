# SpatialPlots_Main 
# the main spatial plotting functions
# ACP 21/4/08
#
# includes 
#   space.plot
#   strata.space.plot
#   scale.plot

`strata.space.plot` <-
function(object,variable,SpaceStrat,TimeStrat,TechStrat,func,nplots=1,multiscale.title,multiscale.cex=1,...)
{
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# strata.space.plot
# function that does stratified spatial plots 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
obj <-missing(object)
vbl <-missing(variable)
spst <-missing(SpaceStrat)
if(any(obj,vbl,spst))stop("function reqiues object variable and SpaceStrat")

if(class(object)=="clData")
{
dataframe <-object@cl
dataname <-"cl"
clnames <-names(object@cl)
if(!(variable %in% clnames))stop("No such variable in the cost object")
}
if(class(object)=="ceData")
{
dataframe <-object@ce
dataname <-"ce"
cenames <-names(object@ce)
if(!(variable %in% cenames))stop("No such variable in the cost object")
}
if(class(object)=="csData")
{
trnames <-names(object@tr)
hhnames <-names(object@hh)
slnames <-names(object@sl)
hlnames <-names(object@hl)
canames <-names(object@ca)
tablenames <-c(rep("tr",length(trnames)),rep("hh",length(hhnames)),
rep("sl",length(slnames)),rep("hl",length(hlnames)),rep("ca",length(canames)))
varindex <-which(is.finite(match(c(trnames,hhnames,slnames,hlnames,canames),variable)))
vartable <-tablenames[varindex]
if(length(vartable)>1) warning("specified variable occurs in more than one csData table")
if(length(vartable)<1) stop("No such variable in the cost object")
newcs <-mergecsData(object)
eval(parse(text=paste("dataframe <-newcs@",vartable[1],sep="")))
dataname <-vartable[1]
if(variable=="lenCls") 
{
plottable <-readline(cat("lenCls occurs in both hl length frequency tables and ca age length tables \n
which do you want to plot? hl or ca \n" ))
eval(parse(text=paste("dataframe <-newcs@",plottable,sep="")))
dataname <-plottable
}
}
if(!(variable %in% names(dataframe)))stop("No such variable in the cost object")
varname <-paste(substitute(variable))
spacename <-paste(substitute(SpaceStrat))
if(!(spacename %in% c("rect","area")))stop("SpaceStat must be either rect or area")

funcname <-paste(substitute(func))

eval(parse(text=paste("variable <-dataframe$",variable,sep="")))
if(is.numeric(variable)==FALSE) stop("You need to spacify a numeric variable")
eval(parse(text=paste("statsqs <-dataframe$",SpaceStrat,sep="")))


default.title <-paste(funcname,"of",varname,"by",spacename,sep=" ")
#print(default.title)
if(nplots>1&missing(multiscale.title))multiscale.title <-default.title


timeperiods <-techtypes <-1
timeindex <-techindex <-1:dim(dataframe)[1]

if(!is.null(TimeStrat))
{
timeperiods <-switch(TimeStrat,year=unique(dataframe$year),quarter=c(1:4),month=c(1:12))
timestrataname <-paste(toupper(substr(TimeStrat,1,1)),substr(TimeStrat,2,nchar(TimeStrat)),sep="")
}



if(!is.null(TechStrat))
{
techtypes <-switch(TechStrat,foCatEu5=unique(dataframe$foCatEu5[!is.na(dataframe$foCatEu5)])
,foCatEu6=unique(dataframe$foCatEu6[!is.na(dataframe$foCatEu6)])
,foCatNat=unique(dataframe$foCatNat[!is.na(dataframe$foCatNat)])
,commCat=unique(dataframe$commCat[!is.na(dataframe$commCat)])
)
}

techtitle <-timetitle <-NULL
# finding the maximum value by strata and setting up a common scale 
maxvals <-NULL
timeindex <-techindex <-1:dim(dataframe)[1]
for(i in 1:length(timeperiods))
{
for(j in 1:length(techtypes))
{
if(!is.null(TimeStrat)) timeindex <-which(dataframe[[TimeStrat]] %in% timeperiods[i])
if(!is.null(TechStrat)) techindex <-which(dataframe[[TechStrat]] %in% techtypes[j])
index <-intersect(timeindex,techindex)
if(length(index)>0) 
{
maxvals <-append(maxvals,as.vector(tapply(variable[index],statsqs[index],func)))
}
}
}
maxstratavalue <-max(maxvals)
commonbreaks <-seq(0,maxstratavalue,length.out=8)

nstrata <-length(timeperiods)*length(techtypes)
pages <-ceiling(nstrata/nplots)
lastpageplot <-nplots*(1:pages)
addtoplot <-c(lastpageplot[lastpageplot<nstrata],nstrata)
addscale <-rep(FALSE,nstrata)
if(!(nplots %in% c(3,5,7,8,10,11)))addscale[addtoplot] <-TRUE


x <-.layout.matrix(nplots)
if(nplots!=1)
{
par(mfrow=c(1,1))
layout(x)
}

# running throught the strat doing a plot for each. 
k <-0
for(i in 1:length(timeperiods))
{
for(j in 1:length(techtypes))
{
k <-k+1
index <-1:dim(dataframe)[1]
plotvars <-space.plot(variable[index],statsqs[index],func,breaks=commonbreaks,plotmap=FALSE,...)
if(!is.null(TimeStrat)) 
{
timeindex <-which(dataframe[[TimeStrat]] %in% timeperiods[i])
timetitle <-paste(timestrataname,"=",ifelse(TimeStrat=="month",month.abb[i],timeperiods[i]),sep=" ")
}
if(!is.null(TechStrat))
{
techindex <-which(dataframe[[TechStrat]] %in% techtypes[j])
techtitle <-paste(TechStrat,"=",techtypes[j],sep=" ")
}
index <-intersect(timeindex,techindex)
if(length(index)>0) #plotmap <-NULL
{
options(warn=-1)
# stops "breaks dont span range of .." warnings when commonbreaks are used
space.plot(variable[index],statsqs[index],func,overlay=T,breaks=commonbreaks,...)
title(sub=paste(timetitle,techtitle,sep=" "))

# adding the multiple scale if its the end of the page
if(k %in% lastpageplot&nplots>1)
{
par(pty="m")
default.mar <-par(mar=c(1,1,1,1))
scale.plot(commonbreaks,place="center",scale.box="n",fillcol=plotvars$cols[1:length(plotvars$cols)],scaletype=plotvars$maptype,multiscale.title,scale.cex=multiscale.cex,cex.max.bubble=plotvars$cex.max.bubble)
scaleplot.mar <-par(mar=default.mar)
}
}
}
}
# resetting graphics paramiters if changed
if(nplots>1)
{
layout(matrix(1,1,1))
par(mfrow=c(1,1))
options(warn=0)
}
}
# ---------------------end of strata.space.plot------------------------------------

 
`space.plot` <-
function(variable,SpaceStrat,func,
xlim=NULL,ylim=NULL,zlim=NULL,xlab=NULL,ylab=NULL,breaks=NULL,
maptype="image",plotmap=TRUE,overlay=FALSE,
ices.divs=FALSE,depths=FALSE,statrects=FALSE,fcoast=FALSE,landmass=FALSE,
pch=1,colour=TRUE,
col.coast="blue",col.cont="grey",col.pch="red",col.rect="grey",col.land="snow2",
col.depth="grey",col.text=1,
scale=FALSE,scale.title="",scale.cex=0.6,scaleplace="bottomright",scale.box="o",
cex.max.bubble=2,threshold=0,digits.text=0,cex.text=1,
...)

{
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# space.plot
# function that takes a numeric variable, a vecotor of 
# statsqs or ICES areas and plots out the 
# spatial information gruoped by func
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(pty="s")
data(NHcoast)
data(finecoast)
data(alldepths)
data(faoAreas)
if(!is.null(breaks)&&substitute(breaks)=="commonbreaks")options(warn=-1)
if(!(maptype %in% c("image","contour","bubble","values"))) 
stop("maptype must be one of: image, contour, bubble, or values")

if(any(is.na(SpaceStrat)))
{
warning("spatial strata included NA's")
newSpaceStrat <-SpaceStrat[!is.na(SpaceStrat)]
variable <-variable[!is.na(SpaceStrat)]
SpaceStrat <-newSpaceStrat
}

if(all(SpaceStrat %in% faoAreas$FAO))
{
#print("areas are FAO")
SpaceStrat <-as.character(faoAreas$ICES[match(SpaceStrat,faoAreas$FAO)])
}

statsqs <-SpaceStrat

if(!is.statsq(SpaceStrat)) statsqs <-convert.icesarea.statsq(SpaceStrat)$statsq
if(!is.statsq(SpaceStrat)&maptype=="contour") stop("Contour plots are not possible for area strata")

# the max plot area
labels <- 0:90
latnames1 <- ifelse(labels < 10, paste("0", labels, sep = ""),as.character(labels))
latnames <-c(rep("-",12),latnames1,rep("-",11))
midlats <-seq(29.5, 86.0, 0.5) + 0.25
lonnames <- c("A0","A1","A2","A3",paste(rep(LETTERS[c(2:8,10:12,13)], rep(10, 11)), rep(0:9,11), sep = ""))
midlons <- (-44:69) + 0.5
 
# The default plotting area covers the limits of the statsq 

latlons <-convert.statsq.lat.lon(statsqs)
navals <-unique(c(which(is.na(latlons$lon)),which(is.na(latlons$lat))))
if(length(navals)>0)
{
latlons$lon <-latlons$lon[-navals]
latlons$lat <-latlons$lat[-navals]
}

#latlons <-latlons[-unique(c(which(is.na(latlons$lon)),which(is.na(latlons$lat)))),]
minx <-min(latlons$lon)
maxx <-max(latlons$lon)
miny <-min(latlons$lat)
maxy <-max(latlons$lat)
if(!is.null(xlim))
{
minx <-xlim[1]
maxx <-xlim[2]
}
if(!is.null(ylim))
{
miny <-ylim[1]
maxy <-ylim[2]
}


xlen <-maxx-minx
ylen <-maxy-miny
midx <-minx+xlen/2
midy <-miny+ylen/2
bigside <-max(xlen,ylen*2)
lowx <-midx-bigside/2
hix <-midx+bigside/2
lowy <-midy-bigside/4
hiy <-midy+bigside/4
areax <-c(lowx-0.5,hix+0.5)
areay <-c(lowy-0.25,hiy+0.25)


# overriding the default x anyd y labels if not specified

if(is.null(xlab)) xlab <-""
if(is.null(ylab)) ylab <-""
if(!(scale.title=="")) scale <-TRUE
if(overlay) par(new=TRUE)
if(!overlay) plot(midlons,midlats,type="n",xlim=areax,ylim=areay,xlab=xlab,ylab=ylab,...)

# calculate statsq values

rects <-expand.grid(lonnames,latnames)
rects1 <-as.factor(paste(rects[,2],rects[,1],sep=""))
valpercell <-tapply(variable,SpaceStrat,func)
statsqs <-as.character(names(valpercell))
values <-as.vector(valpercell)
if(!is.statsq(SpaceStrat))
#if(all(SpaceStrat %in% faoAreas$ICES))
{
casout <-convert.icesarea.statsq(names(valpercell))
statsqs <-casout$statsq
values <-as.vector(valpercell[match(casout$parentarea,toupper(names(valpercell)))])
}
#if(all(SpaceStrat %in% faoAreas$FAO))
#{
#ICESareas <-as.character(faoAreas$ICES[match(names(valpercell),faoAreas$FAO)])
#casout <-convert.icesarea.statsq(ICESareas)
#statsqs <-casout$statsq
#values <-as.vector(valpercell[match(casout$parentarea,toupper(names(valpercell)))])
#}





index <-match(statsqs,rects1)
vals <-rep(NA,length(midlons)^2)
ww <-(1:length(index))[is.finite(index)]
vals[index[ww]] <-(values)[ww]
mat <-matrix(vals,length(midlons),length(midlats))


if (missing(zlim)) 
{
zmin <- min(mat, na.rm = TRUE)
zmax <- max(mat, na.rm = TRUE)
zlim <- c(zmin, zmax)
}
else 
{
mat[mat < zlim[1] & !is.na(mat)] <- zlim[1]
mat[mat > zlim[2] & !is.na(mat)] <- zlim[2]
}

if(!missing(breaks))
{
if(breaks[1]>zlim[1])warning("Breaks dont span the full range of z values")
if(breaks[length(breaks)]<zlim[2])warning("Breaks dont span the full range of z values")
cols <-grey(rev(seq(0.1,0.9,length.out=length(breaks)-1)))
if(colour==T)
{
cols <-rev(heat.colors(length(breaks)-1))
}}

if(missing(breaks))
{
breaks <-seq(zlim[1],zlim[2],length=8)
if(zlim[1]==zlim[2])breaks <-seq(zlim[1],zlim[2],length=2)
cols <-grey(rev(seq(0.1,0.9,length.out=length(breaks)-1)))
if(colour==T)
{
cols <-rev(heat.colors(length(breaks)-1))
}
}

if(!plotmap) scale <-FALSE
# add to the plot

if(maptype=="image"&plotmap)
{
if(all(mat[is.finite(mat)]==0)) cols <-rep("white",length(cols))
image(midlons,midlats,mat,zlim=zlim,add=TRUE,col=cols,breaks=breaks,xlab="",ylab="")

}

if(maptype=="contour"&plotmap)
{
scale <-F
contour(midlons,midlats,mat,col=col.cont,drawlabels=F,nlevels=6,add=TRUE)
}

if(statrects)
{
abline(h=seq(20,90,0.5)-0.5,lty=3,col=col.rect)
abline(v=seq(-50,80,1)+1,lty=3,col=col.rect)
axis(3,midlons,lonnames,cex.axis=0.5,tick=F,line=-0.5)
axis(4,midlats,latnames,cex.axis=0.5,las=3,tick=F,line=-1)

if(maptype=="bubble")
{
if(!fcoast)landmass.polygons("snow2",border=col.coast)
}
}
if(!is.statsq(SpaceStrat)) 
#if(!is.statsq(statsqs)) 
{
ices.divs <-TRUE
landmass.polygons("white",border=col.coast)
}

# adding ices divisions
if(ices.divs)ices.division.lines()
# adding the coastline
if(landmass) fcoast <-FALSE
if(fcoast)lines(finecoast$lon,finecoast$lat,col=col.coast)
if(!(fcoast))lines(NHcoast$lon,NHcoast$lat,col=col.coast)
# adding depths
if(depths[1])
{
for(i in 1:length(depths))
{
ww <-(1:length(alldepths$lon))[depths[i]==alldepths$depth]
lines(alldepths$lon[ww],alldepths$lat[ww],col=col.depth)
}}
if(landmass)landmass.polygons(colour=col.land,border=col.coast)

# the bubble and value plotting comes after coastlines, depths etc so it is not overwriten 
if(maptype=="values"&plotmap)
{
scale <-F
lls <-convert.statsq.lat.lon(as.factor(names(valpercell)))
if(!is.statsq(SpaceStrat)) lls <-convert.icesarea.lat.lon(as.factor(names(valpercell)))
vals <-as.vector(valpercell)
vals <-round(vals,digits.text)
text(lls$lon,lls$lat,vals,col=col.text,cex=cex.text)
}



if(maptype=="bubble"&plotmap)
{
lls <-convert.statsq.lat.lon(as.factor(names(valpercell)))
if(!is.statsq(SpaceStrat)) lls <-convert.icesarea.lat.lon(as.factor(names(valpercell)))
bins <-as.numeric(cut(as.vector(valpercell),breaks,include.lowest=TRUE))
vals <-breaks[bins+1]
size <- cex.max.bubble * ((vals - threshold)/(max(vals,na.rm=T) - threshold))
if(!fcoast)landmass.polygons(col.land,border=col.coast)
points(lls$lon,lls$lat,col=col.pch,cex=size,pch=pch)

}

if(scale)
{
valnames <-breaks
leg <-c(paste(floor(valnames[1:length(valnames)-1]),"<=",floor(valnames[2:(length(valnames))]),sep=" ")) 

if(maptype=="bubble")
{
bins <-as.numeric(cut(as.vector(breaks),breaks,include.lowest=TRUE))
vals <-breaks[bins+1]
scalesize <- cex.max.bubble * ((vals - threshold)/(max(vals,na.rm=T) - threshold))
legend(scaleplace,inset=0.02,leg,pch=rep(pch,length(breaks)),pt.cex=scalesize[2:length(scalesize)]
,cex=scale.cex,title=scale.title,col=col.pch,bty=scale.box)
}


if(maptype=="image")
{
legend(scaleplace,inset=0.02,leg,fill=cols[1:length(cols)],cex=scale.cex,title=scale.title,bty=scale.box)
}
}

box()

out <-list(lon=midlons,lat=midlats,values=mat,bks=breaks,cols=cols,maptype=maptype,scale.cex=scale.cex,cex.max.bubble=cex.max.bubble)
return(invisible(out))

}

#-----------------------end of space.plot----------------------------------------

`scale.plot` <-
function(bks,place="center",multiscale.title="",scaletype,fillcol,scale.box,cex.max.bubble,scale.cex)
{
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function that adds a scale to multiple plots
# produced by strata.space.plot
# needs work as of 21/4/08
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot(1, 1, type='n',xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
valnames <-bks
leg <-c(paste(floor(valnames[1:length(valnames)-1]),"<=",floor(valnames[2:(length(valnames))]),sep=" ")) 
if(scaletype=="image")
{
legend(place,inset=0.02,leg,fill=fillcol,title=multiscale.title,ncol=4,bty=scale.box,cex=scale.cex)
}
if(scaletype=="bubble")
{
threshold <-0
col.pch <-"red"
pch <-1
bins <-as.numeric(cut(as.vector(bks),bks,include.lowest=TRUE))
vals <-bks[bins+1]
scalesize <- cex.max.bubble * ((vals - threshold)/(max(vals,na.rm=T) - threshold))
legend(place,inset=0.02,leg,pch=rep(pch,length(bks)),pt.cex=scalesize[2:length(scalesize)]
,title=multiscale.title,col=col.pch,ncol=4,bty=scale.box,cex=scale.cex)
}
}

#----------------end of scale.plot----------------------------------