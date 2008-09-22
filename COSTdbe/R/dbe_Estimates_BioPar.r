
#########################################################
#  Analytical estimates of empirical weight-at-length,  #
# maturity-at-length, sex-ratio-at-length and variances #
#########################################################



#library(COSTcore)


setGeneric("bpEstim", function(dbeOutput,             #dbeOutput object
                               object,                #csDataCons object
                               ...){

  standardGeneric("bpEstim")

})




setMethod("bpEstim", signature(dbeOutput="dbeOutput",object="csDataCons"), function(dbeOutput,
                                                                                    object,
                                                                                    ...){

As.num <- function(x) as.numeric(as.character(x))
species <- dbeOutput@species
if (all(is.na(species))) stop("no species in 'dbeOutput' object!!")  
if ("all"%in%species) species <- unique(as.character(c(object@ca$spp,object@hl$spp)))

#-------------------------------------------------------------------------------
#restriction of CA and HL tables to specified species (no restriction on catch category, all data is kept)
#-------------------------------------------------------------------------------

ca <- object@ca[object@ca$spp%in%species,]
hl <- object@hl[object@hl$spp%in%species,]
if (nrow(ca)==0) stop("no biological parameters in consolidated object!!")
if (nrow(hl)==0) stop("no hl information!!")


timeStrata <- spaceStrata <- techStrata <- TRUE
if (all(is.na(ca$time))) {timeStrata <- FALSE ; ca$time <- factor("all")}
if (all(is.na(ca$space))) {spaceStrata <- FALSE ; ca$space <- factor("all")}
if (all(is.na(hl$technical))) {techStrata <- FALSE ; hl$technical <- factor("all")}

#-------------------------------------------------------------------------------
#no technical stratification in ca, so, for the homogeneity of dimensions, ca is duplicated within hl$technical field
#-------------------------------------------------------------------------------

CA <- do.call("rbind",lapply( levels(hl$technical),function(x) {df <- ca ; df$technical <- x ; return(df)}))
if (!is.null(CA)) ca <- CA                          

#-------------------------------------------------------------------------------
#recoding of specified biological parameter --> 'bio'
#-------------------------------------------------------------------------------

type <- dbeOutput@param
bio <- switch(type,
              maturity = As.num(cut(As.num(ca$matStage),c(0,1.5,50),labels=c(0,1))),
              sex = As.num(factor(ca$sex,levels=c("F","M"),labels=c("1","0"))),
              weight = As.num(ca$indWt)
              )
if (is.null(bio)) stop("wrong 'param' slot in 'dbeOutput' object!!")
ca$bio <- bio

#-------------------------------------------------------------------------------
#creation of data.frame with numbers of measured fish by length class and strata
#-------------------------------------------------------------------------------

#numbers of measured fish per crossed strata
NL <- tapply(hl$lenNum,list(paste(hl$lenCls,hl$time,hl$space,hl$technical,sep=":-:")),sum,na.rm=TRUE)  #":-:" separator mustn't be in used hl fields
#and we build the data.frame from NL
df <- as.data.frame(t(do.call("cbind",lapply(names(NL),function(x) strsplit(x,":-:")[[1]]))))
df$nL <- NL ; names(df) <- c("lenCls","time","space","technical","nL")

#-------------------------------------------------------------------------------
#creation of data.frame with numbers sub-sampled by length class and strata
#-------------------------------------------------------------------------------

#numbers of measured fish per crossed strata
ML <- tapply(1:nrow(ca),list(paste(ca$lenCls,ca$time,ca$space,ca$technical,sep=":-:")),length)  #":-:" separator mustn't be in used ca fields
#and we build the data.frame from ML
dF <- as.data.frame(t(do.call("cbind",lapply(names(ML),function(x) strsplit(x,":-:")[[1]]))))
dF$mL <- ML ; names(dF) <- c("lenCls","time","space","technical","mL")

#-------------------------------------------------------------------------------
#creation of data.frame with numbers sub-sampled by age and strata
#-------------------------------------------------------------------------------

#numbers of measured fish per crossed strata
MA <- tapply(1:nrow(ca),list(paste(ca$age,ca$time,ca$space,ca$technical,sep=":-:")),length)  #":-:" separator mustn't be in used ca fields
#and we build the data.frame from MA
DF <- as.data.frame(t(do.call("cbind",lapply(names(MA),function(x) strsplit(x,":-:")[[1]]))))
DF$mA <- MA ; names(DF) <- c("age","time","space","technical","mA")

#-------------------------------------------------------------------------------
#df, dF & DF tables are merged to ca
#-------------------------------------------------------------------------------

ca <- merge(ca,df,all.x=TRUE) ; ca$nL[is.na(ca$nL)] <- 0
ca <- merge(ca,dF,all.x=TRUE)
ca <- merge(ca,DF,all.x=TRUE)

#-------------------------------------------------------------------------------
#insertion of r field (r=nL/mL)
#-------------------------------------------------------------------------------

ca$r <- ca$nL/ca$mL

#-------------------------------------------------------------------------------
#sum r of A is put in ca
#-------------------------------------------------------------------------------

rSumA <- tapply(ca$r,list(paste(ca$age,ca$time,ca$space,ca$technical,sep=":-:")),sum,na.rm=TRUE)  #":-:" separator mustn't be in used ca fields
#and we build the data.frame from rSumA
DFA <- as.data.frame(t(do.call("cbind",lapply(names(rSumA),function(x) strsplit(x,":-:")[[1]]))))
DFA$rSumA <- rSumA ; names(DFA) <- c("age","time","space","technical","rSumA")
#and merging is operated
ca <- merge(ca,DFA,all.x=TRUE)

#-------------------------------------------------------------------------------
#calculation of statistical weights w
#-------------------------------------------------------------------------------

ca$w <- ca$mA*ca$r/ca$rSumA




################################################################################
#calculation of mean, var at length, and adjusted mean & var at age (stratified)
################################################################################

MeanAtL <- tapply(ca$bio,list(ca$lenCls,ca$time,ca$space,ca$technical),mean,na.rm=TRUE)
VarAtL <- tapply(ca$bio,list(ca$lenCls,ca$time,ca$space,ca$technical),var,na.rm=TRUE)

MeanAtA <- tapply(ca$w*ca$bio,list(ca$age,ca$time,ca$space,ca$technical),sum,na.rm=TRUE)/tapply(ca$w,list(ca$age,ca$time,ca$space,ca$technical),sum,na.rm=TRUE)
#MeanAtA is inserted in ca for variance calculation
MatA <- tapply(ca$w*ca$bio,list(paste(ca$age,ca$time,ca$space,ca$technical,sep=":-:")),sum,na.rm=TRUE)/tapply(ca$w,list(paste(ca$age,ca$time,ca$space,ca$technical,sep=":-:")),sum,na.rm=TRUE)
DFa <- as.data.frame(t(do.call("cbind",lapply(names(MatA),function(x) strsplit(x,":-:")[[1]]))))
DFa$MatA <- MatA ; names(DFa) <- c("age","time","space","technical","MatA")
#and merging is operated
ca <- merge(ca,DFa,all.x=TRUE)

VarAtA.Num <- tapply(ca$w*((ca$bio-ca$MatA)^2),list(ca$age,ca$time,ca$space,ca$technical),sum,na.rm=TRUE)
VarAtA.Den <- tapply(ca$w,list(ca$age,ca$time,ca$space,ca$technical),sum,na.rm=TRUE)-1  
VarAtA <- VarAtA.Num/VarAtA.Den
#Nan values are in fact 0 values for variance (one fish in the stratum)
VarAtA[is.nan(VarAtA)] <- 0

#dbeOutput must be updated, so result must be formatted
formatt <- function(x,timeStrata,spaceStrata,techStrata,age=FALSE){
column <- expand.grid(dimnames(x))
df <- data.frame(time=column[,2],space=column[,3],technical=column[,4],length=column[,1],value=as.numeric(x))
if (age) names(df)[4] <- "age"
df <- df[!is.na(df$value),] 
if (!timeStrata) df$time <- NA
if (!spaceStrata) df$space <- NA
if (!techStrata) df$technical <- NA
return(df)
}

dbeOutput@lenStruc$estim <- formatt(MeanAtL,timeStrata=timeStrata,spaceStrata=spaceStrata,techStrata=techStrata)
dbeOutput@lenVar <- formatt(VarAtL,timeStrata=timeStrata,spaceStrata=spaceStrata,techStrata=techStrata)
dbeOutput@ageStruc$estim <- formatt(MeanAtA,timeStrata=timeStrata,spaceStrata=spaceStrata,techStrata=techStrata,age=TRUE)
dbeOutput@ageVar <- formatt(VarAtA,timeStrata=timeStrata,spaceStrata=spaceStrata,techStrata=techStrata,age=TRUE)
if (!dbeOutput@methodDesc%in%"analytical") warnings("'methodDesc' slot in 'dbeOutput' object will be updated!!")
dbeOutput@methodDesc <- "analytical"


#-------------------------------------------------------------------------------
#result is returned
#-------------------------------------------------------------------------------

return(dbeOutput)

})






###############
### Example ###
###############

#strDef <- strIni(timeStrata="quarter",spaceStrata="area")
#object <- csDataCons(csDataVal(sole.cs),strDef)
##dbeOutput initial object
#dbeOutput <- dbeObject(desc="My object",species="Solea solea",param="weight",strataDesc=strDef,methodDesc="analytical")
#
#lWeight <- bpEstim(dbeOutput,object)
#dbePlot(lWeight,Slot="lenStruc",step=10,ylab="Mean weight (g)")
#
#dbeOutput2 <- dbeObject(desc="My object",species="Solea solea",param="maturity",strataDesc=strDef,methodDesc="analytical")
#lMaturity <- bpEstim(dbeOutput2,object)
#dbePlot(lMaturity,Slot="lenStruc",step=10,ylab="Maturity ratio")
#
#dbeOutput3 <- dbeObject(desc="My object",species="Solea solea",param="sex",strataDesc=strDef,methodDesc="analytical")
#lSex <- bpEstim(dbeOutput3,object)
#dbePlot(lSex,Slot="lenStruc",step=10,ylab="Sex ratio")
#




