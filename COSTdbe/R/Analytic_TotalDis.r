       #######################
     ###########################
   ###                        ###
 ##### Design-based estimates #####
   ###                        ###
     ###########################
    #############################
    # Tasks related to discards #
    #############################
    


################################################################################
# Internal procedures
################################################################################






#-------------------------------------------------------------------------------
# Formatting function
#-------------------------------------------------------------------------------

procRaise.format <- function(vrbl) {

df <- as.data.frame(do.call("rbind",lapply(names(vrbl),function(x) res <- strsplit(x,":-:")[[1]])))
names(df) <- c("time","space","technical")
df$value <- as.numeric(as.character(vrbl))
return(df)

}







#-------------------------------------------------------------------------------
# Raising methodology based on 4.1 section of 'Raising procedures for discards : Sampling theory' (J. Vigneau, 2006)  --> raising by trip 
#-------------------------------------------------------------------------------

procRaise.trip <- function(csObject,                      #consolidated CS table
                           ceObject,                      #consolidated CE table (same stratification as csObjet)
                           dbeOutp,                       #'dbeOutput' object with descriptive fields 
                           val="weight",                  #value to raise ("weight" or "number")
                           sampPar=TRUE,                  #'sampPar' checks if given species is considered to be automatically sampled
                           ...) {


fraction <- "DIS"

if (dbeOutp@catchCat!=fraction) {
  warnings("'dbeOutput' object doesn't match the method!! 'catchCat' slot will be updated!")
  dbeOutp@catchCat <- "DIS"}
  
if (dbeOutp@methodDesc!="analytical") {
  warnings("'dbeOutput' object doesn't match the method!! 'methodDesc' slot will be updated!")
  dbeOutp@methodDesc <- "analytitcal"}
  
species <- dbeOutp@species
if ("all"%in%species) species <- unique(as.character(csObject@sl$spp))


#according to 'val' and available information, hauls are considered to be sampled or not --> 'sampledFO' method
indSam <- sampledFO(csObject,species=species,fraction=fraction,sampPar=sampPar)
if (val=="weight") indSam <- indSam$sampWt else indSam <- indSam$sampLg
#sampled FOs are ...
samFO <- csObject@hh[!is.na(indSam),c("PSUid","SSUid","time","space","technical")] ; samFO$ind <- 1   #specifies matrix dimensions
#so, values to be considered are ...
if (val=="weight") {
  VAL <- csObject@sl
  VAL$vol <- VAL$wt/1000             # weights in kg
} else {
  VAL <- merge(csObject@hl,csObject@sl[,c("PSUid","SSUid","TSUid","spp","sort","sex","wt","subSampWt")],all.x=TRUE,sort=FALSE)
  VAL$vol <- VAL$lenNum*VAL$wt/VAL$subSampWt
}

VAL <- VAL[VAL$spp%in%species,] ; VAL <- merge(VAL,samFO,all.x=TRUE)
VAL <- VAL[extCatchCat(VAL$sort)%in%fraction & !is.na(VAL$ind),]                #'extCatchCat' function can be found in '0valuesIndex.r' file
if (nrow(VAL)==0) stop("no available sampling data for specified parameters!!")
VALstrat <- paste(VAL$time,VAL$space,VAL$technical,sep=":-:")
samFOstrat <- paste(samFO$time,samFO$space,samFO$technical,sep=":-:")          
HHstrat <- paste(csObject@hh$time,csObject@hh$space,csObject@hh$technical,sep=":-:")

FOind <- tapply(samFO$ind,list(samFO$PSUid,samFO$SSUid,samFOstrat),function(x) return(0))  #index of all sampled FOs

#------
# yij : volume in a haul j of a trip i
#------
if (any(is.na(VAL$vol))) warning("values from hauls defined as sampled hauls are missing!!")
yij <- tapply(VAL$vol,list(factor(VAL$PSUid,levels=levels(factor(samFO$PSUid))),  #all sampled trip should appear
                           factor(VAL$SSUid,levels=levels(factor(samFO$SSUid))),
                           factor(VALstrat,levels=levels(factor(samFOstrat)))),
              sum,na.rm=TRUE)
  #sampled FOs that are not recorded in VAL must also appear in yij as 0 values
yij[is.na(yij) & !is.na(FOind)] <- 0  

#------
# mi : number of sampled hauls in a trip i
#------
mi <- tapply(samFO$ind,list(samFO$PSUid,factor(samFOstrat,levels=dimnames(yij)[[3]])),sum)

#------
# Mi : total number of hauls in a trip i
#------
Mi <- tapply(csObject@hh$SSUid,list(factor(csObject@hh$PSUid,levels=levels(factor(samFO$PSUid))),
             factor(HHstrat,levels=dimnames(yij)[[3]])),function(x) length(unique(x)))

#------
# n : number of sampled trips
#------
n <- tapply(samFO$PSUid,list(factor(samFOstrat,levels=dimnames(yij)[[3]])),function(x) length(unique(x)))

#------
# N : total number of trips
#------
if (all(is.na(ceObject@ce$trpNum))) stop("no available data in 'trpNum' field of 'ceObject' for raising process!!")
if (any(is.na(ceObject@ce$trpNum))) warning("missing values for 'trpNum' field in ceObject!!")
CEstrat <- paste(ceObject@ce$time,ceObject@ce$space,ceObject@ce$technical,sep=":-:")
N <- tapply(ceObject@ce$trpNum,list(factor(CEstrat,levels=dimnames(yij)[[3]])),sum,na.rm=TRUE)

#so, estimate of the total volume is...
yiBar <- apply(yij,c(1,3),sum,na.rm=TRUE)/mi
yI <- N*apply(Mi*yiBar,2,sum,na.rm=TRUE)/n         #if NAs, should be coming from N (unavailable population level data)

#and, its associated variance is...
yiHat <- Mi*yiBar ; yBar <- apply(yiHat,2,sum,na.rm=TRUE)/n
s2i <- apply(aperm(aperm(yij,c(1,3,2))-rep(yiBar,dim(yij)[2]),c(1,3,2))^2,c(1,3),sum,na.rm=TRUE)/(mi-1)
  #Nan & Inf values in 's2i' must be replaced by 0 (mi=1 ==> var=0)
s2i[is.nan(s2i)] <- 0 ; s2i[is.infinite(s2i)] <- 0 
  #first part of the formula of variance
first <- N*N*apply((yiHat-rep(yBar,each=dim(yiHat)[1]))^2,2,sum,na.rm=TRUE)/(n*(n-1))
first[is.nan(first)] <- 0 ; first[is.infinite(first)] <- 0                  #(variance inter trip, so n=1 ==> first=0 ????)
  #second part of the formula
second <- N*apply(Mi*(Mi-mi)*s2i/mi,2,sum,na.rm=TRUE)/n                     #variance intra trip
#VyI
VyI <- first + second


#"dbeOutput" object is updated
est <- switch(val,
        weight="totalW",
        number="totalN")
vr <- switch(val,
        weight="totalWvar",
        number="totalNvar")
        
        
slot(dbeOutp,est)$estim <- procRaise.format(yI)
slot(dbeOutp,vr) <- procRaise.format(VyI)

return(dbeOutp)
}







#-------------------------------------------------------------------------------
# Raising methodology based on 4.2 section of 'Raising procedures for discards : Sampling theory' (J. Vigneau, 2006)  --> raising by hauls 
# >>>----------> Mbar is estimated from the sample (sum(Mi)/n)
#-------------------------------------------------------------------------------

procRaise.fo <- function(csObject,                      #consolidated CS table
                         ceObject,                      #consolidated CE table (same stratification as csObjet)
                         dbeOutp,                       #'dbeOutput' object with descriptive fields 
                         val="weight",                  #value to raise ("weight"or "number")
                         sampPar=TRUE,                  #'sampPar' checks if given species is considered to be automatically sampled
                         ...) {

fraction <- "DIS"

if (dbeOutp@catchCat!=fraction) {
  warnings("'dbeOutput' object doesn't match the method!! 'catchCat' slot will be updated!")
  dbeOutp@catchCat <- "DIS"}
  
if (dbeOutp@methodDesc!="analytical") {
  warnings("'dbeOutput' object doesn't match the method!! 'methodDesc' slot will be updated!")
  dbeOutp@methodDesc <- "analytitcal"}
  
species <- dbeOutp@species
if ("all"%in%species) species <- unique(as.character(csObject@sl$spp))



#according to 'val' and available information, hauls are considered as sampled or not
indSam <- sampledFO(csObject,species=species,fraction=fraction,sampPar=sampPar)
if (val=="weight") indSam <- indSam$sampWt else indSam <- indSam$sampLg
#sampled FOs are ...
samFO <- csObject@hh[!is.na(indSam),c("PSUid","SSUid","time","space","technical")] ; samFO$ind <- 1   #specifies matrix dimensions
#so, values to be considered are ...
if (val=="weight") {
  VAL <- csObject@sl
  VAL$vol <- VAL$wt/1000                    #weights in kg
} else {
  VAL <- merge(csObject@hl,csObject@sl[,c("PSUid","SSUid","TSUid","spp","sort","sex","wt","subSampWt")],all.x=TRUE,sort=FALSE)
  VAL$vol <- VAL$lenNum*VAL$wt/VAL$subSampWt
}

VAL <- VAL[VAL$spp%in%species,] ; VAL <- merge(VAL,samFO,all.x=TRUE)
VAL <- VAL[extCatchCat(VAL$sort)%in%fraction & !is.na(VAL$ind),]
if (nrow(VAL)==0) stop("no available sampling data for specified parameters!!")
#stratification field is built from 'time', 'space' and 'technical' fields for each df
VALstrat <- paste(VAL$time,VAL$space,VAL$technical,sep=":-:")
samFOstrat <- paste(samFO$time,samFO$space,samFO$technical,sep=":-:")          
HHstrat <- paste(csObject@hh$time,csObject@hh$space,csObject@hh$technical,sep=":-:")

FOind <- tapply(samFO$ind,list(samFO$PSUid,samFO$SSUid,samFOstrat),function(x) return(0))  #index of all sampled FOs

#------
# yij : volume in a haul j of a trip i
#------
if (any(is.na(VAL$vol))) warning("values from hauls defined as sampled hauls are missing!!")
yij <- tapply(VAL$vol,list(factor(VAL$PSUid,levels=levels(factor(samFO$PSUid))),  #all sampled trips and FOs should appear
                           factor(VAL$SSUid,levels=levels(factor(samFO$SSUid))),
                           factor(VALstrat,levels=levels(factor(samFOstrat)))),
              sum,na.rm=TRUE)
  #sampled FOs that are not recorded in VAL must also appear in yij as 0 values, so...
yij[is.na(yij) & !is.na(FOind)] <- 0  

#------
# mi : number of sampled hauls in a trip i
#------
mi <- tapply(samFO$ind,list(samFO$PSUid,factor(samFOstrat,levels=dimnames(yij)[[3]])),sum)

#------
# Mi : total number of hauls in a trip i
#------
Mi <- tapply(csObject@hh$SSUid,list(factor(csObject@hh$PSUid,levels=levels(factor(samFO$PSUid))),
             factor(HHstrat,levels=dimnames(yij)[[3]])),function(x) length(unique(x)))

#------
# mBar : mean of number of sampled hauls by trip
#------

mBar <- apply(mi,2,mean,na.rm=TRUE)

#------
# MBar : mean of total number of hauls by trip
#------

MBar <- apply(Mi,2,mean,na.rm=TRUE)

#------
# n : number of sampled trips
#------
n <- tapply(samFO$PSUid,list(factor(samFOstrat,levels=dimnames(yij)[[3]])),function(x) length(unique(x)))

#------
# Mo : total number of FOs
#------
if (all(is.na(ceObject@ce$foNum))) stop("no available data in 'foNum' field of 'ceObject' for raising process!!")
if (any(is.na(ceObject@ce$foNum))) warning("missing values for 'foNum' field in ceObject!!")
CEstrat <- paste(ceObject@ce$time,ceObject@ce$space,ceObject@ce$technical,sep=":-:")
Mo <- tapply(ceObject@ce$foNum,list(factor(CEstrat,levels=dimnames(yij)[[3]])),sum,na.rm=TRUE)     

#so, estimate of the total volume is...
yiBar <- apply(yij,c(1,3),sum,na.rm=TRUE)/mi
yII <- Mo*apply(Mi*yiBar,2,sum,na.rm=TRUE)/apply(Mi,2,sum,na.rm=TRUE)         #if NAs, should be coming from Mo (unavailable population level data)

#and, its associated variance is...
yBarBar <- apply(yiBar,2,sum,na.rm=TRUE)/n                                        
s2i <- apply(aperm(aperm(yij,c(1,3,2))-rep(yiBar,dim(yij)[2]),c(1,3,2))^2,c(1,3),sum,na.rm=TRUE)/(mi-1)
  #Nan & Inf values in 's2i' must be replaced by 0 (mi=1 ==> var=0)
s2i[is.nan(s2i)] <- 0 ; s2i[is.infinite(s2i)] <- 0 
  #first part of the formula
first <- Mo*Mo*apply((yiBar-rep(yBarBar,each=dim(yiBar)[1]))^2,2,sum,na.rm=TRUE)/(n*(n-1))  
first[is.nan(first)] <- 0 ; first[is.infinite(first)] <- 0                  #(variance inter trip, so n=1 ==> first=0 ????)
  #second part of the formula
second <- Mo*Mo*(1-(mBar/MBar))*apply(s2i,2,sum,na.rm=TRUE)/(n*n*mBar)                     #variance intra trip
#VyI
VyII <- first + second

#"dbeOutput" object is updated
est <- switch(val,
        weight="totalW",
        number="totalN")
vr <- switch(val,
        weight="totalWvar",
        number="totalNvar")
        
        
slot(dbeOutp,est)$estim <- procRaise.format(yII)
slot(dbeOutp,vr) <- procRaise.format(VyII)

return(dbeOutp)

}








#-------------------------------------------------------------------------------
# Raising methodology based on 5th section of 'Raising procedures for discards : Sampling theory' (J. Vigneau, 2006)   
# Volume of discards is proportional to an auxiliary variable : Fishing time
#-------------------------------------------------------------------------------

procRaise.time <- function(csObject,                      #consolidated CS table
                           ceObject,                      #consolidated CE table (same stratification as csObjet)
                           dbeOutp,                       #'dbeOutput' object with descriptive fields
                           val="weight",                  #value to raise ("weight" or "number")
                           sampPar=TRUE,                  #'sampPar' checks if given species is considered to be automatically sampled
                           ...) {

fraction <- "DIS"

if (dbeOutp@catchCat!=fraction) {
  warnings("'dbeOutput' object doesn't match the method!! 'catchCat' slot will be updated!")
  dbeOutp@catchCat <- "DIS"}
  
if (dbeOutp@methodDesc!="analytical") {
  warnings("'dbeOutput' object doesn't match the method!! 'methodDesc' slot will be updated!")
  dbeOutp@methodDesc <- "analytitcal"}
  
species <- dbeOutp@species
if ("all"%in%species) species <- unique(as.character(csObject@sl$spp))

#first, 'foDur' information in hh table must be checked
if (any(is.na(csObject@hh$foDur))) warning("Fishing operations with missing 'foDur' information will be considered to be non sampled!!")

#according to 'val' and available information, hauls are considered to be sampled or not
indSam <- sampledFO(csObject,species=species,fraction=fraction,sampPar=sampPar)
if (val=="weight") indSam <- indSam$sampWt else indSam <- indSam$sampLg
#as it's been said in the warning message, hauls with missing fishing duration information are considered to be not sampled
indSam[is.na(csObject@hh$foDur)] <- NA

#sampled FOs are ...
samFO <- csObject@hh[!is.na(indSam),c("PSUid","SSUid","time","space","technical")] ; samFO$ind <- 1   #specifies matrix dimensions
#so, values to be considered are ...
if (val=="weight") {
  VAL <- csObject@sl
  VAL$vol <- VAL$wt/1000          #weights in kg
} else {
  VAL <- merge(csObject@hl,csObject@sl[,c("PSUid","SSUid","TSUid","spp","sort","sex","wt","subSampWt")],all.x=TRUE,sort=FALSE)
  VAL$vol <- VAL$lenNum*VAL$wt/VAL$subSampWt
}


VAL <- VAL[VAL$spp%in%species,] ; VAL <- merge(VAL,samFO,all.x=TRUE)
VAL <- VAL[extCatchCat(VAL$sort)%in%fraction & !is.na(VAL$ind),]
if (nrow(VAL)==0) stop("no available sampling data for specified parameters!!")
#stratification field is built from 'time', 'space' and 'technical' fields for each df
VALstrat <- paste(VAL$time,VAL$space,VAL$technical,sep=":-:")
samFOstrat <- paste(samFO$time,samFO$space,samFO$technical,sep=":-:")          
HHstrat <- paste(csObject@hh$time,csObject@hh$space,csObject@hh$technical,sep=":-:")

FOind <- tapply(samFO$ind,list(samFO$PSUid,samFO$SSUid,samFOstrat),function(x) return(0))  #index of all sampled FOs

#------
# yij : volume in a haul j of a trip i
#------
if (any(is.na(VAL$vol))) warning("values from hauls defined as sampled hauls are missing!!")
yij <- tapply(VAL$vol,list(factor(VAL$PSUid,levels=levels(factor(samFO$PSUid))),  #all sampled trips and FOs should appear
                           factor(VAL$SSUid,levels=levels(factor(samFO$SSUid))),
                           factor(VALstrat,levels=levels(factor(samFOstrat)))),
              sum,na.rm=TRUE)
  #sampled FOs that are not recorded in VAL must also appear in yij as 0 values, so...
yij[is.na(yij) & !is.na(FOind)] <- 0  


#------
# xij : fishing time in a haul j of a trip i
#------
xij <- tapply(csObject@hh$foDur,list(factor(csObject@hh$PSUid,levels=levels(factor(samFO$PSUid))),  
                                     factor(csObject@hh$SSUid,levels=levels(factor(samFO$SSUid))),
                                     factor(HHstrat,levels=levels(factor(samFOstrat)))),
              sum,na.rm=TRUE)
  #sampled FOs are those recorded in yij, so...
xij[is.na(yij)] <- NA 


#------
# mi : number of sampled hauls in a trip i
#------
mi <- tapply(samFO$ind,list(samFO$PSUid,factor(samFOstrat,levels=dimnames(yij)[[3]])),sum)

#------
# Mi : total number of hauls in a trip i
#------
Mi <- tapply(csObject@hh$SSUid,list(factor(csObject@hh$PSUid,levels=levels(factor(samFO$PSUid))),
             factor(HHstrat,levels=dimnames(yij)[[3]])),function(x) length(unique(x)))

#------
# n : number of sampled trips
#------
n <- tapply(samFO$PSUid,list(factor(samFOstrat,levels=dimnames(yij)[[3]])),function(x) length(unique(x)))


#------
# N : total number of trips
#------
if (all(is.na(ceObject@ce$trpNum))) stop("no available data in 'trpNum' field of 'ceObject' for raising process!!")
if (all(is.na(ceObject@ce$foDur))) stop("no available data in 'foDur' field of 'ceObject' for raising process!!")

CEstrat <- paste(ceObject@ce$time,ceObject@ce$space,ceObject@ce$technical,sep=":-:")
if (any(is.na(ceObject@ce$trpNum))) warning("missing values for 'trpNum' field in ceObject!!")
if (any(is.na(ceObject@ce$foDur))) warning("missing values for 'foDur' field in ceObject!!")
N <- tapply(ceObject@ce$trpNum,list(factor(CEstrat,levels=dimnames(yij)[[3]])),sum,na.rm=TRUE)

#------
# X : total stratified fishing time at population level (in minutes)
#------
X <- tapply(ceObject@ce$foDur*60,list(factor(CEstrat,levels=dimnames(yij)[[3]])),sum,na.rm=TRUE)       #in minutes

#so, estimate of the total volume is...
yiBar <- apply(yij,c(1,3),sum,na.rm=TRUE)/mi 
xiBar <- apply(xij,c(1,3),sum,na.rm=TRUE)/mi  
Rhat <- apply(Mi*yiBar,2,sum,na.rm=TRUE)/apply(Mi*xiBar,2,sum,na.rm=TRUE)         
yIII <- X*Rhat                                                                    #if NAs, should be coming from X (unavailable population level data)

#and, its associated variance is...
yiHat <- Mi*yiBar ; xiHat <- Mi*xiBar
s2iPart1 <- yij-xij*rep(Rhat,each=prod(dim(xij)[1:2]))
s2iPart2 <- yiBar-xiBar*rep(Rhat,each=dim(xiBar)[1])
s2i <- apply(aperm(aperm(s2iPart1,c(1,3,2))-rep(s2iPart2,dim(s2iPart1)[2]),c(1,3,2))^2,c(1,3),sum,na.rm=TRUE)/(mi-1)
  #Nan & Inf values in 's2i' must be replaced by 0 (mi=1 ==> var=0)
s2i[is.nan(s2i)] <- 0 ; s2i[is.infinite(s2i)] <- 0 
  #first part of the formula
first <- N*N*apply((yiHat-xiHat*rep(Rhat,each=dim(xiHat)[1]))^2,2,sum,na.rm=TRUE)/(n*(n-1))  
first[is.nan(first)] <- 0 ; first[is.infinite(first)] <- 0                  #not sure about that (variance inter trip, so n=1 ==> first=0 ????)
  #second part of the formula
second <- N*apply(Mi*(Mi-mi)*s2i/mi,2,sum,na.rm=TRUE)/n                     
#VyIII
VyIII <- first + second

#"dbeOutput" object is updated
est <- switch(val,
        weight="totalW",
        number="totalN")
vr <- switch(val,
        weight="totalWvar",
        number="totalNvar")
        
        
slot(dbeOutp,est)$estim <- procRaise.format(yIII)
slot(dbeOutp,vr) <- procRaise.format(VyIII)

return(dbeOutp)

}







#-------------------------------------------------------------------------------
# Raising methodology based on 5th section of 'Raising procedures for discards : Sampling theory' (J. Vigneau, 2006) 
# Volume of discards is proportional to an auxiliary variable : Landings of the species (species to be specified with 'landSpp' parameter)
#-------------------------------------------------------------------------------

procRaise.landings <- function(csObject,                      #consolidated CS table
                               ceObject,                      #consolidated CE table (same stratification as csObjet)
                               clObject,                      #consolidated CL table (same stratification as csObjet)
                               dbeOutp,                       #'dbeOutput' object with descriptive fields
                               landSpp=as.character(NA),      #NA (ie landSpp=species),"all" (raising variable is complete landings), or a character vector with species name(s)
                               val="weight",                  #value to raise ("weight"or "number")
                               sampPar=TRUE,                  #'sampPar' checks if given species is considered as automatically sampled
                               ...) {


if (dbeOutp@catchCat!="DIS") {
  warnings("'dbeOutput' object doesn't match the method!! 'catchCat' slot will be updated!")
  dbeOutp@catchCat <- "DIS"}
  
if (dbeOutp@methodDesc!="analytical") {
  warnings("'dbeOutput' object doesn't match the method!! 'methodDesc' slot will be updated!")
  dbeOutp@methodDesc <- "analytitcal"}
  
species <- dbeOutp@species
if ("all"%in%species) species <- unique(as.character(csObject@sl$spp))

if (all(is.na(landSpp))) landSpp <- species
if ("all"%in%landSpp) landSpp <- species

#according to 'val' and available information, hauls are considered to be sampled or not.
#Since both fraction must be sampled in that case,...
indSamDis <- sampledFO(csObject,species=species,fraction="DIS",sampPar=sampPar)
if (val=="weight") indSamDis <- indSamDis$sampWt else indSamDis <- indSamDis$sampLg
indSamLan <- sampledFO(csObject,species=landSpp,fraction="LAN",sampPar=sampPar)
if (val=="weight") indSamLan <- indSamLan$sampWt else indSamLan <- indSamLan$sampLg
indSam <- (!is.na(indSamDis) & !is.na(indSamLan))   #only hauls where discards and landings have been sampled are kept for the calculation 

#sampled FOs are ...
samFO <- csObject@hh[indSam,c("PSUid","SSUid","time","space","technical")] ; samFO$ind <- 1   #specifies matrix dimensions
#so, values to be considered are ...
if (val=="weight") {
  VAL <- csObject@sl
  VAL$vol <- VAL$wt/1000         #weights in kg
} else {
  VAL <- merge(csObject@hl,csObject@sl[,c("PSUid","SSUid","TSUid","spp","sort","sex","wt","subSampWt")],all.x=TRUE,sort=FALSE)
  VAL$vol <- VAL$lenNum*VAL$wt/VAL$subSampWt
}

LAN <- VAL
VAL <- VAL[VAL$spp%in%species,] ; VAL <- merge(VAL,samFO,all.x=TRUE)
VAL <- VAL[extCatchCat(VAL$sort)%in%"DIS" & !is.na(VAL$ind),]
if (nrow(VAL)==0) stop("no available discards sampling data for specified parameters!!")
#stratification field is built from 'time', 'space' and 'technical' fields for each df
VALstrat <- paste(VAL$time,VAL$space,VAL$technical,sep=":-:")
samFOstrat <- paste(samFO$time,samFO$space,samFO$technical,sep=":-:")          
HHstrat <- paste(csObject@hh$time,csObject@hh$space,csObject@hh$technical,sep=":-:")

#sampled landed weights for specified species
LAN <- LAN[LAN$spp%in%landSpp,] ; LAN <- merge(LAN,samFO,all.x=TRUE)
LAN <- LAN[extCatchCat(LAN$sort)%in%"LAN" & !is.na(LAN$ind),]
if (nrow(LAN)==0) stop("no available landings sampling data for specified parameters!!")
LANstrat <- paste(LAN$time,LAN$space,LAN$technical,sep=":-:")



#index of all sampled FOs
FOind <- tapply(samFO$ind,list(samFO$PSUid,samFO$SSUid,samFOstrat),function(x) return(0))  

#------
# yij : volume in a haul j of a trip i
#------
if (any(is.na(VAL$vol))) warning("values from hauls defined as sampled hauls are missing!!")
yij <- tapply(VAL$vol,list(factor(VAL$PSUid,levels=levels(factor(samFO$PSUid))),  #all sampled trips and FOs should appear
                           factor(VAL$SSUid,levels=levels(factor(samFO$SSUid))),
                           factor(VALstrat,levels=levels(factor(samFOstrat)))),
              sum,na.rm=TRUE)
  #sampled FOs that are not recorded in VAL must also appear in yij as 0 values, so...
yij[is.na(yij) & !is.na(FOind)] <- 0  


#------
# xij : landed weight of specified species in a haul j of a trip i  (in kg)
#------
xij <- tapply(LAN$wt/1000,list(factor(LAN$PSUid,levels=levels(factor(samFO$PSUid))),  
                               factor(LAN$SSUid,levels=levels(factor(samFO$SSUid))),
                               factor(LANstrat,levels=levels(factor(samFOstrat)))),
              sum,na.rm=TRUE)
  #sampled FOs are those recorded in yij, so...
xij[is.na(yij)] <- NA 


#------
# mi : number of sampled hauls in a trip i
#------
mi <- tapply(samFO$ind,list(samFO$PSUid,factor(samFOstrat,levels=dimnames(yij)[[3]])),sum)

#------
# Mi : total number of hauls in a trip i
#------
Mi <- tapply(csObject@hh$SSUid,list(factor(csObject@hh$PSUid,levels=levels(factor(samFO$PSUid))),
             factor(HHstrat,levels=dimnames(yij)[[3]])),function(x) length(unique(x)))

#------
# n : number of sampled trips
#------
n <- tapply(samFO$PSUid,list(factor(samFOstrat,levels=dimnames(yij)[[3]])),function(x) length(unique(x)))


#------
# N : total number of trips
#------
CEstrat <- paste(ceObject@ce$time,ceObject@ce$space,ceObject@ce$technical,sep=":-:")
if (any(is.na(ceObject@ce$trpNum))) warning("missing values for 'trpNum' field in ceObject!!")
N <- tapply(ceObject@ce$trpNum,list(factor(CEstrat,levels=dimnames(yij)[[3]])),sum,na.rm=TRUE)

#------
# X : total stratified landed weights at population level for landSpp species (in kg)
#------
CLt <- clObject@cl[clObject@cl$taxon%in%landSpp,]
CLstrat <- paste(CLt$time,CLt$space,CLt$technical,sep=":-:")
if (any(is.na(CLt$landWt))) warning("missing values for 'landWt' field in clObject!!")
X <- tapply(CLt$landWt,list(factor(CLstrat,levels=dimnames(yij)[[3]])),sum,na.rm=TRUE)     

#so, estimate of the total volume is...
yiBar <- apply(yij,c(1,3),sum,na.rm=TRUE)/mi 
xiBar <- apply(xij,c(1,3),sum,na.rm=TRUE)/mi  
Rhat <- apply(Mi*yiBar,2,sum,na.rm=TRUE)/apply(Mi*xiBar,2,sum,na.rm=TRUE)         #if NaN, no landings and no discards  --> 0 values
Rhat[is.nan(Rhat)] <- 0
yIII <- X*Rhat                                                                    #if NAs, should be coming from X (unavailable population level data)

#and, its associated variance is...
yiHat <- Mi*yiBar ; xiHat <- Mi*xiBar
s2iPart1 <- yij-xij*rep(Rhat,each=prod(dim(xij)[1:2]))
s2iPart2 <- yiBar-xiBar*rep(Rhat,each=dim(xiBar)[1])
s2i <- apply(aperm(aperm(s2iPart1,c(1,3,2))-rep(s2iPart2,dim(s2iPart1)[2]),c(1,3,2))^2,c(1,3),sum,na.rm=TRUE)/(mi-1)
  #Nan & Inf values in 's2i' must be replaced by 0 (mi=1 ==> var=0)
s2i[is.nan(s2i)] <- 0 ; s2i[is.infinite(s2i)] <- 0 
  #first part of the formula
first <- N*N*apply((yiHat-xiHat*rep(Rhat,each=dim(xiHat)[1]))^2,2,sum,na.rm=TRUE)/(n*(n-1))  
first[is.nan(first)] <- 0 ; first[is.infinite(first)] <- 0                  #variance inter trip, so n=1 ==> first=0 ????
  #second part of the formula
second <- N*apply(Mi*(Mi-mi)*s2i/mi,2,sum,na.rm=TRUE)/n                     
#VyIII
VyIII <- first + second
          
#"dbeOutput" object is updated
est <- switch(val,
        weight="totalW",
        number="totalN")
vr <- switch(val,
        weight="totalWvar",
        number="totalNvar")
        
        
slot(dbeOutp,est)$estim <- procRaise.format(yIII)
slot(dbeOutp,vr) <- procRaise.format(VyIII)

return(dbeOutp)
}








#-------------------------------------------------------------------------------
# New raising methodology based on fishing days at population level
#-------------------------------------------------------------------------------

procRaise.fd <- function(csObject,                      #consolidated CS table
                         ceObject,                      #consolidated CE table (same stratification as csObjet)
                         dbeOutp,                       #'dbeOutput' object with descriptive fields
                         val="weight",                  #value to raise ("weight"or "number")
                         sampPar=TRUE,                  #'sampPar' checks if given species is considered as automatically sampled
                         ...) {
    
fraction <- "DIS"

if (dbeOutp@catchCat!=fraction) {
  warnings("'dbeOutput' object doesn't match the method!! 'catchCat' slot will be updated!")
  dbeOutp@catchCat <- "DIS"}
  
if (dbeOutp@methodDesc!="analytical") {
  warnings("'dbeOutput' object doesn't match the method!! 'methodDesc' slot will be updated!")
  dbeOutp@methodDesc <- "analytitcal"}
  
species <- dbeOutp@species
if ("all"%in%species) species <- unique(as.character(csObject@sl$spp))
                  
#according to 'val' and available information, hauls are considered to be sampled or not
indSam <- sampledFO(csObject,species=species,fraction=fraction,sampPar=sampPar)
if (val=="weight") indSam <- indSam$sampWt else indSam <- indSam$sampLg
#sampled FOs are ...
samFO <- csObject@hh[!is.na(indSam),c("PSUid","SSUid","time","space","technical","date")] ; samFO$ind <- 1   #specifies matrix dimensions
#unused levels are deleted
samFO <- model.frame(PSUid~SSUid+time+space+technical+date+ind,data=samFO,drop.unused.level=TRUE)
#so, values to be considered are ...
if (val=="weight") {
  VAL <- csObject@sl
  VAL$vol <- VAL$wt/1000             # weights in kg
} else {
  VAL <- merge(csObject@hl,csObject@sl[,c("PSUid","SSUid","TSUid","spp","sort","sex","wt","subSampWt")],all.x=TRUE,sort=FALSE)
  VAL$vol <- VAL$lenNum*VAL$wt/VAL$subSampWt
}

VAL <- VAL[VAL$spp%in%species,] ; VAL <- merge(VAL,samFO,all.x=TRUE)
VAL <- VAL[extCatchCat(VAL$sort)%in%fraction & !is.na(VAL$ind),]
if (nrow(VAL)==0) stop("no available sampling data for specified parameters!!")
VALstrat <- paste(VAL$time,VAL$space,VAL$technical,sep=":-:")
samFOstrat <- paste(samFO$time,samFO$space,samFO$technical,sep=":-:")          
HHstrat <- paste(csObject@hh$time,csObject@hh$space,csObject@hh$technical,sep=":-:")

FOind <- tapply(samFO$ind,list(samFO$PSUid,samFO$date,samFO$SSUid,samFOstrat),function(x) return(0))  #index of all sampled FOs

#------
# yikj : volume in a haul j of a fishing day k, in a trip i
#------
if (any(is.na(VAL$vol))) warning("values from hauls defined as sampled hauls are missing!!")
yikj <- tapply(VAL$vol,list(factor(VAL$PSUid,levels=levels(factor(samFO$PSUid))),  #all sampled trip should appear
                           factor(VAL$date,levels=levels(factor(samFO$date))),
                           factor(VAL$SSUid,levels=levels(factor(samFO$SSUid))),
                           factor(VALstrat,levels=levels(factor(samFOstrat)))),
              sum,na.rm=TRUE)
  #sampled FOs that are not recorded in VAL must also appear in yikj as 0 values
yikj[is.na(yikj) & !is.na(FOind)] <- 0  

#------
# mik : number of sampled hauls in a fishing day k of a trip i
#------
mik <- tapply(samFO$ind,list(samFO$PSUid,samFO$date,factor(samFOstrat,levels=dimnames(yikj)[[4]])),sum)

#------
# Mik : total number of hauls in a fishing day k of a trip i
#------
Mik <- tapply(csObject@hh$SSUid,list(factor(csObject@hh$PSUid,levels=levels(factor(samFO$PSUid))),
             factor(csObject@hh$date,levels=levels(factor(samFO$date))),
             factor(HHstrat,levels=dimnames(yikj)[[4]])),function(x) length(unique(x)))

#------
# di : number of sampled fishing day in a trip i
#------
di <- tapply(samFO$date,list(samFO$PSUid,factor(samFOstrat,levels=dimnames(yikj)[[4]])),function(x) length(unique(x)))

#------
# Di : total number of fishing day in a trip i
#------
Di <- tapply(csObject@hh$date,list(factor(csObject@hh$PSUid,levels=levels(factor(samFO$PSUid))),
             factor(HHstrat,levels=dimnames(yikj)[[4]])),function(x) length(unique(x)))

#------
# n : number of sampled trips
#------
n <- tapply(samFO$PSUid,list(factor(samFOstrat,levels=dimnames(yikj)[[4]])),function(x) length(unique(x)))

#------
# D : total number of fishing days
#------
CEstrat <- paste(ceObject@ce$time,ceObject@ce$space,ceObject@ce$technical,sep=":-:")
if (any(is.na(ceObject@ce$daysAtSea))) warning("missing values for 'daysAtSea' field in ceObject!!")
D <- tapply(ceObject@ce$daysAtSea,list(factor(CEstrat,levels=dimnames(yikj)[[4]])),sum,na.rm=TRUE)

#so, estimate of the total volume is...
yikBar <- apply(yikj,c(1,2,4),sum,na.rm=TRUE)/mik
yikHat <- Mik*yikBar
yiBar <- apply(yikHat,c(1,3),sum,na.rm=TRUE)/di
yBarBar <- apply(Di*yiBar,2,sum,na.rm=TRUE)/apply(Di,2,sum,na.rm=TRUE)
yIV <- D*apply(Di*yiBar,2,sum,na.rm=TRUE)/apply(Di,2,sum,na.rm=TRUE)         #if NAs, should be coming from D (unavailable population level data)

#and, its associated variance is...
  #first part of the formula of variance
first <- n*((D/apply(Di,2,sum,na.rm=TRUE))^2)*(1-apply(Di,2,sum,na.rm=TRUE)/D)*apply(aperm(aperm(yiBar,c(2,1))-yBarBar,c(2,1))^2,2,sum,na.rm=TRUE)/(n-1)
  #Nan & Inf values in 'first' must be replaced by 0 (mi=1 ==> var=0)
first[is.nan(first)] <- 0 ; first[is.infinite(first)] <- 0 

  #second part of the formula of variance
second <- D*apply(Di*(Di-di)*apply(aperm(aperm(yikHat,c(1,3,2))-as.numeric(yiBar),c(1,3,2))^2,c(1,3),sum,na.rm=TRUE)/(di*(di-1)),2,sum,na.rm=TRUE)/apply(Di,2,sum,na.rm=TRUE)
  #Nan & Inf values in 'second' must be replaced by 0 (mi=1 ==> var=0)
second[is.nan(second)] <- 0 ; second[is.infinite(second)] <- 0 

  #third part of the formula of variance
stepp <- Mik*(Mik-mik)*apply(aperm(aperm(yikj,c(1,2,4,3))-as.numeric(yikBar),c(1,2,4,3))^2,c(1,2,4),sum,na.rm=TRUE)/(mik*(mik-1))
third <- D*apply(Di*apply(stepp,c(1,3),sum,na.rm=TRUE)/di,2,sum,na.rm=TRUE)/apply(Di,2,sum,na.rm=TRUE)
  #Nan & Inf values in 'third' must be replaced by 0 (mi=1 ==> var=0)
third[is.nan(third)] <- 0 ; third[is.infinite(third)] <- 0 

  
#VyIV
VyIV <- first + second + third

#"dbeOutput" object is updated
est <- switch(val,
        weight="totalW",
        number="totalN")
vr <- switch(val,
        weight="totalWvar",
        number="totalNvar")
        
        
slot(dbeOutp,est)$estim <- procRaise.format(yIV)
slot(dbeOutp,vr) <- procRaise.format(VyIV)

return(dbeOutp)                        
}



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################










################################################################################
# Methods to be exported
################################################################################


setGeneric("totVolume", function(dbeOutput,
                                 csObject,
                                 ceObject,
                                 clObject,
                                 ...){
	standardGeneric("totVolume")}
)



setMethod("totVolume", signature(dbeOutput="dbeOutput",csObject="csDataCons",ceObject="ceDataCons",clObject="missing"), function(dbeOutput,
                                                                                                                                 csObject,
                                                                                                                                 ceObject,
                                                                                                                                 type="trip",   #or "fo", "fd", "time"
                                                                                                                                 val="weight",  #or "number"
                                                                                                                                 sampPar=TRUE,
                                                                                                                                 ...){
if (type=="landings") stop("'landings' type requires a cs, a ce and a cl object!!")
eval(parse('',text=paste("procRaise.",type,"(csObject,ceObject,dbeOutput,val=val,sampPar=sampPar)",sep="")))                                                                                                           

})


setMethod("totVolume", signature(dbeOutput="dbeOutput",csObject="csDataCons",ceObject="ceDataCons",clObject="clDataCons"), function(dbeOutput,
                                                                                                                                    csObject,
                                                                                                                                    ceObject,
                                                                                                                                    clObject,
                                                                                                                                    landSpp=as.character(NA),
                                                                                                                                    val="weight",  #or "number"
                                                                                                                                    sampPar=TRUE,
                                                                                                                                    ...){

procRaise.landings(csObject,ceObject,clObject,dbeOutput,landSpp=landSpp,val=val,sampPar=sampPar)                                                                                                          

})










###############
###############
##  Example  ##
###############
###############


#library(COSTcore)
#
#source("C:/Documents and Settings/mmerzere/Bureau/0valuesIndex.r")
#
##consolidated datasets are built for testing 
#strDef <- strIni(timeStrata="quarter",techStrata="foCatEu5")
#csObject <- csDataCons(csDataVal(sole.cs),strDef)
#clObject <- clDataCons(clDataVal(sole.cl),strDef)
#ceObject <- ceDataCons(ceDataVal(sole.ce),strDef)
#
##random foNum are created in ceObject
#ceObject@ce$foNum <- sample(c(15:30),nrow(ceObject@ce),replace=TRUE)*ceObject@ce$trpNum           
#
##random foDur are created in ceObject
#ceObject@ce$foDur <- sample(c(60:180),nrow(ceObject@ce),replace=TRUE)*ceObject@ce$trpNum            
#
##random daysAtSea values are created in ceObject
#ceObject@ce$daysAtSea <- sample(c(1:5),nrow(ceObject@ce),replace=TRUE)*ceObject@ce$trpNum           
#
############################################
#
##dbeOutput initial object
#obj <- dbeObject(desc="My object",species="Solea solea",catchCat="DIS",strataDesc=strDef,methodDesc="analytical")
#obj
#
##raising by trip
#newObj1 <- totVolume(obj,csObject,ceObject)
##raising by haul
#newObj2 <- totVolume(obj,csObject,ceObject,type="fo")
##raising by fishing day
#newObj3 <- totVolume(obj,csObject,ceObject,type="fd")
##raising by fishing duration
#newObj4 <- totVolume(obj,csObject,ceObject,type="time")
##raising by total landings
#newObj5 <- totVolume(obj,csObject,ceObject,clObject)
#



