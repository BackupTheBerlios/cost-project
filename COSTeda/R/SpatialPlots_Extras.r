# SpatialPlots
# functions for spatial plotting 
# ACP 21/4/08
# includes:  
#   mergecsData
#   layout.matrix
#   landmass.polygons
#   is.statsq
#   convert.statsq.lat.lon
#   convert.icesarea.statsq
#   convert.icesarea.lat.lon
#   convert.statsq.icesarea
#   ices.division.lines
#   ices.division.names




`mergecsData` <-
function(csobj)
{
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function that merges the tables of a csData object
# so that $rect, $area, $date, $foCatNat, 
# foCatEu5 and foCatEu6 are appended to each of the 
# tables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(class(csobj)%in%c("csData","csDataVal")!=TRUE) stop("this function only works on a csData object")


ca6colstring <-apply(csobj@ca[,c(1:6)],1,paste,collapse=".")
tr6colstring <-apply(csobj@tr[,c(1:6)],1,paste,collapse=".")
hh6colstring <-apply(csobj@hh[,c(1:6)],1,paste,collapse=".")
hh7colstring <-apply(csobj@hh[,c(1:7)],1,paste,collapse=".")
sl7colstring <-apply(csobj@sl[,c(1:7)],1,paste,collapse=".")
hl7colstring <-apply(csobj@hl[,c(1:7)],1,paste,collapse=".")


hlindex <-match(hl7colstring,hh7colstring)
if(any(is.na(hlindex))) warning("The key fields between hh and hl tables dont all match")
csobj@hl$rect <-csobj@hh$rect[hlindex]
csobj@hl$area <-csobj@hh$area[hlindex]
csobj@hl$date <-csobj@hh$date[hlindex]
csobj@hl$foCatNat <-csobj@hh$foCatNat[hlindex]
csobj@hl$foCatEu5 <-csobj@hh$foCatEu5[hlindex]
csobj@hl$foCatEu6 <-csobj@hh$foCatEu6[hlindex]
csobj@hl$quarter <-as.numeric(substr(quarters(as.POSIXlt(csobj@hl$date)),2,3))
csobj@hl$month <-as.POSIXlt(csobj@hl$date)$mon+1
csobj@hl$yearfromdate <-as.numeric(substr(csobj@hl$date,1,4))



slindex <-match(sl7colstring,hh7colstring)
if(any(is.na(slindex))) warning("The key fields between hh and sl tables dont all match")
csobj@sl$rect <-csobj@hh$rect[slindex]
csobj@sl$area <-csobj@hh$area[slindex]
csobj@sl$date <-csobj@hh$date[slindex]
csobj@sl$foCatNat <-csobj@hh$foCatNat[slindex]
csobj@sl$foCatEu5 <-csobj@hh$foCatEu5[slindex]
csobj@sl$foCatEu6 <-csobj@hh$foCatEu6[slindex]
csobj@sl$quarter <-as.numeric(substr(quarters(as.POSIXlt(csobj@sl$date)),2,3))
csobj@sl$month <-as.POSIXlt(csobj@sl$date)$mon+1
csobj@sl$yearfromdate <-as.numeric(substr(csobj@sl$date,1,4))


trindex <-match(tr6colstring,hh6colstring)
if(any(is.na(trindex))) warning("The key fields between tr and hh tables dont all match")
csobj@tr$rect <-csobj@hh$rect[trindex]
csobj@tr$area <-csobj@hh$area[trindex]
csobj@tr$date <-csobj@hh$date[trindex]
csobj@tr$foCatNat <-csobj@hh$foCatNat[trindex]
csobj@tr$foCatEu5 <-csobj@hh$foCatEu5[trindex]
csobj@tr$foCatEu6 <-csobj@hh$foCatEu6[trindex]
csobj@tr$quarter <-as.numeric(substr(quarters(as.POSIXlt(csobj@tr$date)),2,3))
csobj@tr$month <-as.POSIXlt(csobj@tr$date)$mon+1
csobj@tr$yearfromdate <-as.numeric(substr(csobj@tr$date,1,4))
 
caindex <-match(ca6colstring,tr6colstring)
if(any(is.na(caindex))) warning("The key fields between tr and ca tables dont all match")
#csobj@ca$rect <-csobj@tr$rect[caindex]
#csobj@ca$area <-csobj@tr$area[caindex]
#csobj@ca$date <-csobj@tr$date[caindex]
csobj@ca$foCatNat <-csobj@tr$foCatNat[caindex]
csobj@ca$foCatEu5 <-csobj@tr$foCatEu5[caindex]
csobj@ca$foCatEu6 <-csobj@tr$foCatEu6[caindex]
#csobj@ca$quarter <-as.numeric(substr(quarters(as.POSIXlt(csobj@ca$date)),2,3))
#csobj@ca$yearfromdate <-as.numeric(substr(csobj@ca$date,1,4))
#csobj@ca$monthfromdate <-as.numeric(substr(csobj@ca$date,6,7))

#csobj@hh$quarter <-as.numeric(substr(quarters(as.POSIXlt(csobj@hh$date)),2,3))
#csobj@hh$month <-as.POSIXlt(csobj@hh$date)$mon+1
#csobj@hh$yearfromdate <-as.numeric(substr(csobj@hh$date,1,4))

return(csobj)
}
#------------end of mergecsData-------------------------------

`.layout.matrix` <-
function(nplots)
{
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# layout.matrix
# function that sets the parameters for layout depending 
# on the number of plots required per page. 
# Adds one to nplots for the multiple scale to be plotted 
# along the bottom of the page used in strata.space.plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(!(nplots %in% c(1,2,4,6,9,12)))stop ("nplots must be one of 1,2,4,6,9,or 12")
if(nplots==4)
{
x <-matrix(0,7,6)
x[1:3,1:3] <-1
x[1:3,4:6] <-2
x[4:6,1:3] <-3
x[4:6,4:6] <-4
x[7,1:6] <-5
}

if(nplots==2)
{
x <-matrix(0,4,6)
x[1:3,1:3] <-1
x[1:3,4:6] <-2
x[4,1:6] <-3
}

if(nplots==1)
{
x <-matrix(1,1,1)
}

if(nplots==6)
{
x <-matrix(0,5,6)
x[1:2,1:2] <-1
x[1:2,3:4] <-2
x[1:2,5:6] <-3
x[3:4,1:2] <-4
x[3:4,3:4] <-5
x[3:4,5:6] <-6
x[5,1:6] <-7
}

if(nplots==9)
{
x <-matrix(0,7,6)
x[1:2,1:2] <-1
x[1:2,3:4] <-2
x[1:2,5:6] <-3
x[3:4,1:2] <-4
x[3:4,3:4] <-5
x[3:4,5:6] <-6
x[5:6,1:2] <-7
x[5:6,3:4] <-8
x[5:6,5:6] <-9
x[7,1:6] <-10
}

if(nplots==12)
{
x <-matrix(0,7,8)
x[1:2,1:2] <-1
x[1:2,3:4] <-2
x[1:2,5:6] <-3
x[1:2,7:8] <-4
x[3:4,1:2] <-5
x[3:4,3:4] <-6
x[3:4,5:6] <-7
x[3:4,7:8] <-8
x[5:6,1:2] <-9
x[5:6,3:4] <-10
x[5:6,5:6] <-11
x[5:6,7:8] <-12
x[7,1:8] <-13
}
return(x)
}

#---------------------end of layout matrix---------------------------

`landmass.polygons` <-
function(colour="lightgrey",border="grey",...)
{
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# landmass.polygons
# function that adds landmass polygons
# for the NOAA 1:5000000 coastline using the coordinates
# defined by landmass in the list "landmasses"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data(NHcoast)
data(landmasses)
NHcoast$string[24171:24335] <-722

alllandmasses <-c("africa","andros","anglesey","arran","baliaricsC","baliaricsW","baltic1",
     "baltic2","baltic3","baltic4","bear","cephalonia","chios","corsica",
"crete","cyprus","dutch1","dutch2","dutch3","dutch4",
"eire","euboea","europe",
"far1","far2","far3","far4","far5","far6",
"greenland","greenland1","greenland2","greenland3","greenland4"
,"greenland5","greenland6","greenland7","greenland8","greenland9"
,"greenland10","greenland11","greenland12","greenland13","greenland14","greenland15",
"greenland16","greenland17",
"iceland","islay","jura","lemnos","lesvos",
"lewis","malta","man","mull","naxos","nor1","nor2","nor3","nor4","nor5","nor6",
"nor7","nor8","nor9","nor10","nor11","nor12","nor13",
"nor14","nor15","nor16","nor17","nor18","nor19","nor20","nor21",
"novaya","novaya1","novaya2","ork","ork2",
"peloponnesus","rhodes","sardinia","shet","sicily","skye",
"spit1","spit2","spit3","spit4","spit5","spit6","spit7","uist","ukmainland"
,"white") 
       
for(i in 1:length(alllandmasses))
{
eval(parse(text=paste("lon <-landmasses$",alllandmasses[i],"$lon",sep="")))
eval(parse(text=paste("lat <-landmasses$",alllandmasses[i],"$lat",sep="")))
polygon(lon,lat,col=colour,border=border,...)
}
box()
}

#------------------enf of landmass polygons-------------------------------

`is.statsq` <-
function(x)
{
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# is.statsq
# returns TRUE 
# if x is 4 characters long, 3rd of which in a capital letter 
# and the other 3 values are 0 to 9 characters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(class(x)!="character") x <-as.character(x)
numberchars <-as.character(0:9)
all(c(all(nchar(x)==4),
all(substr(x,3,3) %in% LETTERS),
all((substr(x,4,4) %in% numberchars)),
all((substr(x,1,1) %in% numberchars)),
all((substr(x,2,2) %in% numberchars))
))
}
#--------------end of is.statsq-------------------------

`convert.statsq.lat.lon` <-
function (statsq) 
{
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# convert.statsq.lat.lon
# function that returns the centres of statsqs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    part1 <- substr(statsq, 1, 2)
    part2 <- substr(statsq, 3, 4)
    labels <- 0:90
    latlabels <- ifelse(labels < 10, paste("0", labels, sep = ""), 
        as.character(labels))
    latvalues <- seq(35.5, 80.5, 0.5) + 0.25
   lonlabels <- c("A0","A1","A2","A3",paste(rep(LETTERS[c(2:8,10:12,13)], rep(10, 11)), rep(0:9, 
        11), sep = ""))
lonvalues <- (-44:69) + 0.5
    indx <- match(part1, latlabels)
    lat <- latvalues[indx]
    indx <- match(part2, lonlabels)
    lon <- lonvalues[indx]
if (any(is.na(lat)) | any(is.na(lon))) 
warning("Some stat squares have not been recognised.")
out <-list(lon=lon,lat=lat)
    return(out)
}
#-----------------end of convert.statsq.lat.lon----------------

`convert.icesarea.statsq` <-
function(areas)
{
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# convert.icesarea.statsq
# function that returns the constituent statsqs of an  
# of an ICES area 
# needs data frame ICESAreaRects
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data(ICESAreaRects)
areas <-toupper(as.character(areas))
uniarea <-unique(areas)
a <-ICESAreaRects
a$statsq <-toupper(as.character(a$statsq))
a$division <-toupper(as.character(a$division))
a$subdivision <-toupper(as.character(a$subdivision))
a$subarea <-toupper(as.character(a$subarea))
b <-stack(a,select=c(subarea,division,subdivision))
b$statsq <-rep(a$statsq,3)
statsqs <-as.character(b$statsq[which(!is.na(match(b$values,areas)))])
parentarea <-as.character(b$values[which(!is.na(match(b$values,areas)))])
out <-list(statsq=statsqs,parentarea=parentarea)
return(out)
}

#-----------------end of convert.icesarea.statsq------------------------

`convert.icesarea.lat.lon` <-
function(icesarea)
{
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function that returns the approx centre 
# of an ICES area
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lons <-c(3,3,2.5,-9,-15,-15,-15,-10.5,-19,3,-5,-11,-10.5,-7,-7,-3.5,0,3,-18,-35,-28,-28,-3.5,
-2,-8,-6,-14.25,-14.25,-10,40)
lats <-c(52,55.5,59.5,58,58,53.5,50.5,61.5,67,67,53.75,53.5,50.5,51,49,49.5,50.25,78,72.5,62,54,43,47
,45,46,44,46,40,40,78)
areas <-c("IVc","IVb","IVa","VIa","VIb","VIIc","VIIk","Vb","Va","IIa","VIIa","VIIb","VIIj","VIIg",
"VIIh","VIIc","VIId","IIb","XIVa","XIVb","XII","X","VIIIa","VIIIb","VIIId","VIIIc","VIIIe"
,"IXb","IXa","I")

index <-match(icesarea,areas)
if(any(is.na(index)))warning("some of the ices areas not recognised")
out <-list(lon=lons[index],lat=lats[index])
return(out)
}

#-----------------end of convert.icesarea.lat.lon------------------------

convert.statsq.icesarea <-function(statsqs)
{
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# convert.statsq.area
# function that gives the ICES subarea, division and subdivision
# for a given statsq
# needs ICESAreaRects
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data(ICESAreaRects)
aa <-ICESAreaRects
if(is.statsq(statsqs)!=TRUE)stop("statistical rectangles not in expected format")
check <-match(statsqs,aa$statsq)
if(any(is.na(check)!=FALSE))warning("some of your statsqs are not recognised")
subdivs <-as.character(aa$subdivision[match(statsqs,aa$statsq)])
divs <-as.character(aa$division[match(statsqs,aa$statsq)])
subarea <-as.character(aa$subarea[match(statsqs,aa$statsq)])
out <-list(subdivs=subdivs,divs=divs,subarea=subarea)
return(out)
}


#-------------------end of convert.statsq.icesarera
`ices.division.lines` <-
function (division = NULL, area = NULL, lty = 1, col = 1, lwd = 1) 
{
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ices.division.lines
# function that draws ICES division and subarea boundaries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    division.list <- c("I","IIa", "IIb", "IIIa","IIIb","IIIc", "IVa", "IVb", "IVc", 
        "Va", "Vb", "VIa", "VIb", "VIIa", "VIIb", "VIIc", "VIId", 
        "VIIe", "VIIf", "VIIg", "VIIh", "VIIj", "VIIk", "VIIIa", 
        "VIIIb", "VIIIc", "VIIId", "VIIIe", "IXa", "IXb", "X", 
        "XII", "XIVa", "XIVb")
    area.list <- c("I","II", "II", "III", "IV", "IV", "IV", "V", 
        "V", "VI", "VI", "VII", "VII", "VII", "VII", "VII", "VII", 
        "VII", "VII", "VII", "VII", "VII", "VIII", "VIII", "VIII", 
        "VIII", "IX", "IX", "X", "XII", "XIV", "XIV")
    if (is.null(division) & is.null(area)) 
        division <- division.list
    if (!is.null(area)) {
        div <- division.list[area.list %in% area]
        if (is.null(division)) 
            division <- div
        else division <- c(division, div)
    }
    missing.division <- division[!(division %in% division.list)]
    if (length(missing.division) > 0) 
        warning(paste("we don't cater for these divisions:", 
            missing.division))
  if("I" %in% division)
{
lat <-c(90,90,77.5,NA,71,72,72,90)
lon <-c(30,69,69,NA,28,28,30,30)
lines(lon, lat, lty = lty, col = col, lwd = lwd)
}  




if ("IIa" %in% division) {
        lat <- c(73.5, 73.5, 72, 72, 71, NA, 62, 62, 63, 63, 
            73.5)
        lon <- c(-11, 30, 30, 28, 28, NA, 5.2, -4, -4, -11, -11)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("IIb" %in% division) {
        lat <- c(90, 90, 73.5, 73.5)
        lon <- c(-11, 30, 30, -11)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("IIIa" %in% division) {
        lat <- c(58.2, 57.5, 57.5, 57, 57, NA)
        lon <- c(7, 7, 8, 8, 8.2, NA)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("IVa" %in% division & !("IVb" %in% division)) {
        lat <- c(62, 62, NA, 58.2, 57.5, 57.5, NA, 58.5, 62)
        lon <- c(-4, 5.2, NA, 7, 7, -1.7, NA, -4, -4)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("IVb" %in% division & !("IVa" %in% division)) {
        lat <- c(57.5, 57.5, 57, 57, NA, 53.5, 53.5)
        lon <- c(-1.7, 8, 8, 8.5, NA, 0, 7)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("IVa" %in% division & "IVb" %in% division) {
        lat <- c(62, 62, NA, 58, 57.5, NA, 58.5, 62, NA, 57.5, 
            57.5, 57, 57, NA, 53.5, 53.5)
        lon <- c(-4, 5, NA, 7, 7, NA, -4, -4, NA, -1.7, 8, 8, 
            8.5, NA, 0, 7)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("IVc" %in% division) {
        lat <- c(53.5, 53.5, NA, 51, 51)
        lon <- c(0, 7, NA, 1, 2)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("Va" %in% division) {
        lat <- c(68, 68, 63, 63, 62, 62)
        lon <- c(-27, -11, -11, -15, -15, -27)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("Vb" %in% division) {
        lat <- c(60, 60, 60.5, 60.5, 63, 63, 60)
        lon <- c(-15, -5, -5, -4, -4, -15, -15)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("VIa" %in% division) {
        lat <- c(60, 60, 60.5, 60.5, 58.5, NA, 55, 55, NA, 54.5, 
            54.5, 60)
        lon <- c(-12, -5, -5, -4, -4, NA, -5, -6, NA, -8, -12, 
            -12)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("VIb" %in% division) {
        lat <- c(60, 60, 54.5, 54.5, 60)
        lon <- c(-18, -12, -12, -18, -18)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("VIIa" %in% division) {
        lat <- c(55, 55, NA, 52, 52)
        lon <- c(-6, -5.2, NA, -4.7, -7.5)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("VIIb" %in% division) {
        lat <- c(54.5, 54.5, NA, 52.5, 52.5)
        lon <- c(-12, -8, NA, -9.5, -12)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("VIIc" %in% division) {
        lat <- c(54.5, 54.5, 52.5, 52.5, 54.5)
        lon <- c(-18, -12, -12, -18, -18)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("VIId" %in% division) {
        lat <- c(51, 51, NA, 49.5, 50.7)
        lon <- c(1, 2, NA, -2, -2)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("VIIe" %in% division) {
        lat <- c(49.5, 50.7, NA, 50, 50, 49.5, 49.5, 48, 48)
        lon <- c(-2, -2, NA, -5.5, -7, -7, -5, -5, -4.5)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("VIIf" %in% division) {
        lat <- c(NA, 51.5, 51, 51, 50.5, 50.5, 50, 50)
        lon <- c(NA, -5, -5, -6, -6, -7, -7, -5.5)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("VIIg" %in% division) {
        lat <- c(52, 52, NA, 51.5, 51, 51, 50.5, 50.5, 50, 50, 
            51.5)
        lon <- c(-7.5, -4.7, NA, -5, -5, -6, -6, -7, -7, -9, 
            -9)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("VIIh" %in% division) {
        lat <- c(50, 50, 49.5, 49.5, 48, 48, 50)
        lon <- c(-9, -7, -7, -5, -5, -9, -9)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("VIIj" %in% division) {
        lat <- c(52.5, 52.5, NA, 51.5, 48, 48, 50)
        lon <- c(-12, -9.5, NA, -9, -9, -12, -12)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("VIIk" %in% division) {
        lat <- c(52.5, 52.5, 48, 48, 52.5)
        lon <- c(-18, -12, -12, -18, -18)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("VIIIa" %in% division) {
        lat <- c(48, 48, 47.5, 47.5, 47, 47, 46, 46)
        lon <- c(-4.5, -8, -8, -6, -6, -5, -5, -0.9)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("VIIIb" %in% division) {
        lat <- c(46, 46, 45.5, 45.5, 44.5, 44.5, 43.5)
        lon <- c(-0.9, -4, -4, -3, -3, -2, -2)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("VIIIc" %in% division) {
        lat <- c(44.5, 44.5, 43.3, NA, 43, 43, 43, 48)
        lon <- c(-11, -2, -2, NA, -9, -11, -18, -18)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("VIIId" %in% division) {
        lat <- c(48, 48, 47.5, 47.5, 47, 47, 46, 46, 45.5, 45.5, 
            44.5, 44.5, 48)
        lon <- c(-11, -8, -8, -6, -6, -5, -5, -4, -4, -3, -3, 
            -11, -11)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("VIIIe" %in% division) {
        lat <- c(48, 48, 43, 43, 48)
        lon <- c(-18, -11, -11, -18, -18)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("IXa" %in% division) {
        lat <- c(43, 43, NA, 36, 36, 43)
        lon <- c(-11, -9, NA, -5.4, -11, -11)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("IXb" %in% division) {
        lat <- c(43, 43, 36, 36, 43)
        lon <- c(-18, -11, -11, -18, -18)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("X" %in% division) {
        lat <- c(48, 48, 36, 36,48)
        lon <- c(-42, -18, -18, -42,-42)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("XII" %in% division) {
        lat <- c(62, 62, 60, 60, 48, 48, 59, 59, 62)
        lon <- c(-27, -15, -15, -18, -18, -42, -42, -27, -27)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("XIVa" %in% division) {
        lat <- c(83.5, 90, 90, 68, 68, 68.7)
        lon <- c(-40, -40, -11, -11, -27, -27)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
    }
    if ("XIVb" %in% division) {
        lat <- c(68.7, 59, 59, 60)
        lon <- c(-27, -27, -44, -44)
        lines(lon, lat, lty = lty, col = col, lwd = lwd)
        }
    if("IIIa" %in% division) {
     
lon <-c(10.74488,11.34143)
lat <-c(56.15972,55.93602)
lines(lon,lat,lty = lty, col = col, lwd = lwd)
lines(c(12.43510,12.23625),c(56.20944,56.08516),lty = lty, col = col, lwd = lwd)
         }
if("IIIb" %in% division) {    
lines(c(12.43510,12.78308),c(55.31461,55.36433),lty = lty, col = col, lwd = lwd)        
lines(c(12.43510,12.23625),c(56.20944,56.08516),lty = lty, col = col, lwd = lwd)        
        
        }
if("IIIc" %in% division) {    
lon <-c(10.74488,11.34143)
lat <-c(56.15972,55.93602)
lines(lon,lat,lty = lty, col = col, lwd = lwd)
lines(c(11.88826,11.88826),c(54.54407,54.17123),lty = lty, col = col, lwd = lwd)
 }
    
}



#--------------------------end of ices.division.lines-------------------

ices.division.names <-function(text.cex=1)
{
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ices.division.names
# function that adds ICES area names to a plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(cex=text.cex)
text(c(3,3,2.5,-9,-15,-15,-15,-10.5,-19,3,-5,-11,-10.5,-7,-7,-3.5,0,3,-18,-35,-28,-28,-3.5,
-2,-8,-6,-14.25,-14.25,-10,40),
c(52,55.5,59.5,58,58,53.5,50.5,61.5,67,67,53.75,53.5,50.5,51,49,49.5,50.25,78,72.5,62,54,43,47
,45,46,44,46,40,40,78),
c("IVc","IVb","IVa","VIa","VIb","VIIc","VIIk","Vb","Va","IIa","VIIa","VIIb","VIIj","VIIg",
"VIIh","VIIc","VIId","IIb","XIVa","XIVb","XII","X","VIIIa","VIIIb","VIIId","VIIIc","VIIIe"
,"IXb","IXa","I"))
par(cex=1)
}

#-------------end of ices.division.names----------------------------------


