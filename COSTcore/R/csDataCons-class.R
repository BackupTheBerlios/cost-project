#====================================================================
#
# EJ, 24/09/2007
# csData-class
#
# Data hierarchy
# tr --> hh --> sl --> hl
# tr --> ca
#  
#====================================================================

#====================================================================
# Abbreviates
#
# samp = sampling
# country = ctry
# flag =flg 
# vessel = vsl
# trip = trp
# number = num
# project = proj
# power = pwr
# fishing operation = fo
# method = meth
# depth = dep
# selectivity =sel
# device =dev
# commercial = comm
# weight = wt
# class = cls
# otolith = oto
#====================================================================

#====================================================================
# Class definition and validity check
#====================================================================

valcscData <- function(object){
	tr <- object@tr
	hh <- object@hh
	sl <- object@sl
	hl <- object@hl
	ca <- object@ca

	# I will rely o the prototype to check col names and size. I'm not sure it's a good strategy !
	obj <- new("csDataCons")
	tr0 <- obj@tr
	hh0 <- obj@hh
	sl0 <- obj@sl
	hl0 <- obj@hl
	ca0 <- obj@ca
	
	# check columns
	if(checkNms(tr, names(tr0))==FALSE) stop("Check slot candidate \"tr\" columns' size and names.")
	if(checkNms(hh, names(hh0))==FALSE) stop("Check slot candidate \"hh\" columns' size and names.")
	if(checkNms(sl, names(sl0))==FALSE) stop("Check slot candidate \"sl\" columns' size and names.")
	if(checkNms(hl, names(hl0))==FALSE) stop("Check slot candidate \"hl\" columns' size and names.")
	if(checkNms(ca, names(ca0))==FALSE) stop("Check slot candidate \"ca\" columns' size and names.")
	
	# check PK (ToDo)

	# check data integrity
#	if(checkDataIntegrity(tr[,1:6], hh[,1:6])==FALSE) stop("Data integrity problem in table \"hh\". Missing related records in \"tr\".")
#	if(checkDataIntegrity(hh[,1:7], sl[,1:7])==FALSE) stop("Data integrity problem in table \"sl\". Missing related records in \"hh\".")
#	if(checkDataIntegrity(sl[,1:14], hl[,1:14])==FALSE) stop("Data integrity problem in table \"hl\". Missing related records in \"sl\".")
#	if(checkDataIntegrity(tr[,1:6], ca[,1:6])==FALSE) stop("Data integrity problem in table \"ca\". Missing related records in \"tr\".")

	# Everything is fine
	return(TRUE)
}

setClass("csDataCons",
	representation(
		desc="character",
		tr="data.frame",
		hh="data.frame",
		sl="data.frame",
		hl="data.frame",
		ca="data.frame"
	),
	prototype(
		desc="my stock",
		tr=data.frame(
			PSUid=as.numeric(NA), # field to be created by the EDA process, one per each record = PK
			time=as.factor(NA),
			space=as.factor(NA),
			technical=as.factor(NA),
			sampType=as.factor(NA), # PK
			landCtry=as.factor(NA), # PK
			vslFlgCtry=as.factor(NA), # PK
			#year=as.numeric(NA), # PK 	=> time
			proj=as.factor(NA), # PK
			trpCode=as.factor(NA), # PK
			#vslLen=as.numeric(NA), 	=> tech
			#vslPwr=as.numeric(NA), 	=> tech 
			#vslSize=as.numeric(NA), 	=> tech 
			#vsType=as.character(NA), 	=> tech
			foNum=as.numeric(NA), 
			daysAtSea=as.numeric(NA), 
			vslId=as.numeric(NA), 
			sampCtry=as.factor(NA), 
			sampMeth=as.factor(NA)),
		hh=data.frame(
			PSUid=as.numeric(NA), # FK
			SSUid=as.numeric(NA), # field to be created by the EDA process, one per each record = FK+PK
			time=as.factor(NA),
			space=as.factor(NA),
			technical=as.factor(NA),
			sampType=as.character(NA), # FK
			landCtry=as.character(NA), # FK
			vslFlgCtry=as.character(NA), # FK
			#year=as.numeric(NA), # PK 	=> time
			proj=as.character(NA), # FK
			trpCode=as.character(NA), # FK
			staNum=as.numeric(NA), # PK
			foVal=as.character(NA),
			aggLev=as.character(NA),
			catReg=as.character(NA),
			sppReg=as.character(NA),
			#date=as.character(NA), 	=> time
			#time=as.character(NA), 	=> time
			foDur=as.numeric(NA),
			latIni=as.numeric(NA),
			lonIni=as.numeric(NA),
			latFin=as.numeric(NA),
			lonFin=as.numeric(NA),
			#area=as.character(NA),	=> space
			#rect=as.character(NA),	=> space
			foDep=as.numeric(NA)
			#waterDep=as.numeric(NA), => space
			#foCatNat=as.character(NA),	=> tech
			#foCatEu5=as.character(NA),	=> tech
			#foCatEu6=as.character(NA),	=> tech
			#gear=as.character(NA),	=> tech
			#meshSize=as.numeric(NA),	=> tech
			#selDev=as.character(NA),	=> tech
			#meshSizeSelDev=as.numeric(NA)	=> tech
			),
		sl=data.frame(
			PSUid=as.numeric(NA), # FK
			SSUid=as.numeric(NA), # field to be created by the EDA process, one per each record = FK+PK
			TSUid=as.numeric(NA), # field to be created by the EDA process, one per each record = FK+PK
			time=as.factor(NA),
			space=as.factor(NA),
			technical=as.factor(NA),
			sort=as.factor(NA),
			sampType=as.factor(NA), # FK
			landCtry=as.factor(NA), # FK
			vslFlgCtry=as.factor(NA), # FK
			#year=as.numeric(NA), # PK 	=> time
			proj=as.factor(NA), # FK
			trpCode=as.factor(NA), # FK
			staNum=as.numeric(NA), # FK
			spp=as.factor(NA), # PK 
			#catchCat=as.character(NA), # PK	=> sort 
			#landCat=as.character(NA), # PK	=> sort 
			#commCatScl=as.character(NA), # PK	=> sort
			#commCat=as.character(NA), # PK	=> sort
			#subSampCat=as.character(NA), # PK	=> sort
#			valCode=as.factor(NA), 
			wt=as.numeric(NA), 
			subSampWt=as.numeric(NA), 
			lenCode=as.factor(NA)),
		hl=data.frame(
			PSUid=as.numeric(NA), # This field helps on linking with TR table
			SSUid=as.numeric(NA), # This field helps on linking with HH table 
			TSUid=as.numeric(NA), # FK
			time=as.factor(NA),
			space=as.factor(NA),
			technical=as.factor(NA),
			sort=as.factor(NA),
			sampType=as.character(NA), # FK
			landCtry=as.character(NA), # FK
			vslFlgCtry=as.character(NA), # FK
			#year=as.numeric(NA), # PK 	=> time
			proj=as.character(NA), # FK
			trpCode=as.character(NA), # FK
			staNum=as.numeric(NA), # FK
			spp=as.character(NA), # FK 
			sex=as.character(NA), # PK
			#catchCat=as.character(NA), # PK	=> sort 
			#landCat=as.character(NA), # PK	=> sort 
			#commCatScl=as.character(NA), # PK	=> sort
			#commCat=as.character(NA), # PK	=> sort
			#subSampCat=as.character(NA), # PK	=> sort
			lenCls=as.numeric(NA), # PK
			lenNum=as.numeric(NA)),
		ca=data.frame(
			PSUid=as.numeric(NA), # FK
			SSUid=as.numeric(NA), # must match SSUid in HH so that info about tech can be used if necessary
			time=as.factor(NA),
			space=as.factor(NA),
			technical=as.factor(NA),
			sort=as.factor(NA),
			sampType=as.character(NA), # FK
			landCtry=as.character(NA), # FK
			vslFlgCtry=as.character(NA), # FK
			#year=as.numeric(NA), # FK	=> time
			#quarter=as.numeric(NA),	=> time
			#month=as.numeric(NA),	=> time
			proj=as.character(NA), # FK
			trpCode=as.character(NA), # FK
			staNum=as.numeric(NA), # PK
			spp=as.character(NA), # PK 
			sex=as.character(NA), # PK
			#catchCat=as.character(NA), # PK	=> sort 
			#landCat=as.character(NA), # PK	=> sort 
			#commCatScl=as.character(NA), # PK	=> sort
			#commCat=as.character(NA), # PK	=> sort
			stock=as.character(NA), # PK
			#area=as.character(NA),	=> space
			#rect=as.character(NA),	=> space
			lenCls=as.numeric(NA), # PK
			age=as.numeric(NA), # PK
			fishId=as.numeric(NA), # PK
			lenCode=as.character(NA),
			plusGrp=as.character(NA),
			otoWt=as.numeric(NA),
			otoSide=as.character(NA),
			indWt=as.numeric(NA),
			matScale=as.character(NA),
			matStage=as.character(NA))
	),
	validity=valcscData
)

##====================================================================
## Class constructor
##====================================================================
#setGeneric("csDataCons", function(object, ...){
#	standardGeneric("csDataCons")
#	}
#)
#
#setMethod("csDataCons", signature("csDataVal"), function(object, desc="Unknown stock", ...){
#
#	#------------------------------------------------------------------------------
#	# create SUid
#	#------------------------------------------------------------------------------
#	suid <- createSUid(object)
#
#	#------------------------------------------------------------------------------
#	# time
#	#------------------------------------------------------------------------------
#	# hh
#	hh <- suid$hh
#	tm <- as.Date(hh$date)
#	tm <- paste(hh$year,quarters(tm),sep=".")
#	hh$time <- tm
#	# tr
#	tm <- unique(suid$hh[,c("PSUid", "time")])
#	tr <- merge(suid$tr, tm, all.x=TRUE)
#	# sl
#	sl <- suid$sl
#	sl$time <- hh$time[match(sl$SSUid,hh$SSUid)]
#	# hl
#	hl <- suid$hl
#	hl$time <- hh$time[match(hl$SSUid,hh$SSUid)]
#	# ca
#	ca <- suid$ca
#	ca$time <- paste(ca$year, paste("Q", ca$quarter, sep=""), sep=".")
#	
#	#------------------------------------------------------------------------------
#	# tech
#	#------------------------------------------------------------------------------
#	te1 <- apply(tr[,c("vslLen", "vslPwr", "vslSize", "vslType")], 1,paste, collapse=".")
#	te2 <- apply(hh[,c("foCatNat","foCatEu5","foCatEu6","gear","meshSize","selDev","meshSizeSelDev")], 1,paste, collapse=".")
#	#hh
#	hh$technical <- paste(te1[match(hh$PSUid,tr$PSUid)], te2, sep="-")
#	# tr
#	te <- unique(hh[,c("PSUid", "technical")])
#	tr <- merge(tr, te, all.x=TRUE)
#	# sl
#	sl$technical <- hh$technical[match(sl$SSUid,hh$SSUid)]
#	# hl
#	hl$technical <- hh$technical[match(hl$SSUid,hh$SSUid)]
#	# ca
#	ca$technical <- hh$technical[match(ca$SSUid,hh$SSUid)]
#	
#	# Not relevant
#	
#	#------------------------------------------------------------------------------
#	# space
#	#------------------------------------------------------------------------------
#	# hh
#	hh$space <- apply(hh[,c("area","rect","waterDep")], 1,paste, collapse=".") 
#	# tr
#	# not relevant
#	# sl
#	sl$space <- hh$space[match(sl$SSUid,hh$SSUid)]
#	# hl
#	hl$space <- hh$space[match(hl$SSUid,hh$SSUid)]
#	# ca
#	ca$space <- apply(ca[,c("area","rect")], 1,paste, collapse=".")
#	
#	#------------------------------------------------------------------------------
#	# sort
#	#------------------------------------------------------------------------------
#	# sl
#	sl$sort <- apply(sl[,c("catchCat","landCat","commCatScl", "commCat", "subSampCat")], 1,paste, collapse=".")
#	# hl
#	hl$sort <- sl$sort[match(hl$TSUid,sl$TSUid)]
#	# ca
#	ca$sort <- apply(ca[,c("catchCat","landCat","commCatScl", "commCat")], 1,paste, collapse=".")
#	
#	#------------------------------------------------------------------------------
#	# create csDataCons
#	#------------------------------------------------------------------------------
#	csc <- csDataCons()
#	tr <- tr[,match(names(tr(csc)),names(tr))]
#	hh <- hh[,match(names(hh(csc)),names(hh))]
#	sl <- sl[,match(names(sl(csc)),names(sl))]
#	hl <- hl[,match(names(hl(csc)),names(hl))]
#	ca <- ca[,match(names(ca(csc)),names(ca))]
#	new("csDataCons", tr=tr, hh=hh, sl=sl, hl=hl, ca=ca)
#})
#
#
#setMethod("csDataCons", signature("missing"), function(desc="Unknown stock", ...){
#	# create object and name columns properly 
#	new("csDataCons", desc=desc)
#})
#

setGeneric("csDataCons", function(object,objStrat,...){
	standardGeneric("csDataCons")
})


setMethod("csDataCons", signature("csDataVal","strIni"), function(object,
                                                                  objStrat,
                                                                  desc="Unknown stock",  
                                                                  ...){  

timeStrata <- objStrat@timeStrata 
spaceStrata <- objStrat@spaceStrata 
techStrata <- objStrat@techStrata 
#sorting <- objStrat@sorting
tpRec <- objStrat@tpRec
spRec <- objStrat@spRec
tcRec <- objStrat@tcRec


TR <- object@tr 
HH <- object@hh 
SL <- object@sl
HL <- object@hl 
CA <- object@ca

#-------------------------------------------------------------------------------
# Special case if tech strata = commercial cat
# this inforation stored in sl must be written in hh$technical
# and as many rows in hh as there are "commCat" modalities in sl must be created
#-------------------------------------------------------------------------------
cc <- is.na(techStrata)==FALSE & techStrata=="commCat"
if (cc) {
  HH$ind <- 1:nrow(HH)
  HH <- merge(HH,unique(SL[,c("sampType",
                              "landCtry",
                              "vslFlgCtry",
                              "year",
                              "proj",
                              "trpCode",
                              "staNum",
                              "commCat")]),sort=FALSE,all=TRUE)
  HH <- HH[order(HH$ind),]}


#-------------------------------------------------------------------------------
# Creation of the 3 stratification fields in hh
#-------------------------------------------------------------------------------

# Creation of a data.frame containing the different modalities of time stratification
month = sapply(HH$date,function(x) as.numeric(strsplit(x,"-")[[1]][2]))
time.DF <- data.frame(month = month, quarter = ceiling(month/3), semester = ceiling(month/6),year=HH$year)


if (is.na(timeStrata)) {
  HH$time <- NA 
  tpRec <- as.list(NA)
} else {
  HH$time <- time.DF[,timeStrata]}  
     
if (is.na(spaceStrata)) {
  HH$space <- NA 
  spRec <- as.list(NA)
} else {
  HH$space <- HH[,spaceStrata]}
  
if (is.na(techStrata)) {
  HH$technical <- NA 
  tcRec <- as.list(NA)
} else {
  HH$technical <- HH[,techStrata]}


#-------------------------------------------------------------------------------
# Recoding the 3 new stratification fields following user post-stratification          <<<- there's surely a more simple way to do this
#-------------------------------------------------------------------------------
if (!is.na(tpRec[1])) {
  Typ <- class(HH$time) 
  hhT <- factor(HH$time) 
  Lev <- levels(hhT)[!levels(hhT)%in%tpRec$from]
  HH$time <- factor(hhT,levels=c(Lev,tpRec$from),labels=c(Lev,tpRec$to))
  eval(parse('',text=paste("HH$time <- as.",Typ,"(as.character(HH$time))",sep="")))}
  
if (!is.na(spRec[1])) {
  Typ <- class(HH$space) 
  hhS <- factor(HH$space) 
  Lev <- levels(hhS)[!levels(hhS)%in%spRec$from]
  HH$space <- factor(hhS,levels=c(Lev,spRec$from),labels=c(Lev,spRec$to)) 
  eval(parse('',text=paste("HH$space <- as.",Typ,"(as.character(HH$space))",sep="")))}
  
if (!is.na(tcRec[1])) {
  Typ <- class(HH$technical) 
  hhC <- factor(HH$technical)
  Lev <- levels(hhC)[!levels(hhC)%in%tcRec$from]
  HH$technical <- factor(hhC,levels=c(Lev,tcRec$from),labels=c(Lev,tcRec$to))
  eval(parse('',text=paste("HH$technical <- as.",Typ,"(as.character(HH$technical))",sep="")))}
                        

#-------------------------------------------------------------------------------
# Creation of the fields PSUid and SSUid in hh
#-------------------------------------------------------------------------------
#PSUid is defined by the combination of a trip, time, space and technical strata
psuid <- apply(HH[,c("sampType","landCtry","vslFlgCtry","proj","trpCode","time","space","technical")],1,paste,collapse=":-:")     #delim modified
HH$PSUid <- factor(psuid,levels=unique(psuid),labels=1:length(unique(psuid)))
HH$PSUid <- as.numeric(as.character(HH$PSUid))
#SSUid
HH <- HH[order(HH$PSUid),]                                    #Need to be ordered like the result of SSUid below                                         
HH$SSUid <- unlist(tapply(HH$PSUid,list(HH$PSUid),function(x) 1:length(x)))



#-------------------------------------------------------------------------------
# Update of the fields PSUid, time, space and technical in tr from hh
# and recalculation of number of days at sea and number of FO following the new stratification
#-------------------------------------------------------------------------------
tr <- merge(TR,unique(HH[,c("sampType",
                            "landCtry",
                            "vslFlgCtry",
                            "year",
                            "proj",
                            "trpCode",
                            "PSUid",
                            "time",
                            "space",
                            "technical")]),sort=FALSE,all=TRUE)
#'number of days at sea' & 'number of FO' fields are only to be updated for trips with aggLev="H" so... 
HHh <- HH[HH$aggLev=="H",]
#table of updated 'number of days' and 'number of FO' fields 
count <- function(x){length(unique(x))}
df <- data.frame(
   PSUid=tapply(HHh$PSUid,as.character(HHh$PSUid),unique),
   nOP=tapply(HHh$staNum,as.character(HHh$PSUid),count),
   nDay=tapply(HHh$date,as.character(HHh$PSUid),count)
   )
#substitution of 'number of days' and 'number of FO' fields by updated values in tr table 
index <- match(tr$PSUid,df$PSUid)
index <- index[!is.na(index)]
tr$foNum[index] <- df$nOP
tr$daysAtSea[index] <- df$nDay


                        
#-------------------------------------------------------------------------------
# Addition of the fields PSUid, SSUid, time, space and technical in sl from hh
# and creation of TSUid and sorting fields
#-------------------------------------------------------------------------------
#if tech strata = commercial cat, technical field must be included at first
if (cc) SL$technical <- SL$commCat

sl <- merge(SL,HH[,c("sampType",
                     "landCtry",
                     "vslFlgCtry",
                     "year",
                     "proj",
                     "trpCode",
                     "staNum",
                     "time",
                     "space",
                     "technical",
                     "PSUid",
                     "SSUid")],sort=FALSE,all.x=TRUE)

# creation of 'sort' field 
fields <- c("catchCat","landCat","commCatScl","commCat","subSampCat")
sl$sort <- apply(sl[,fields],1,paste,collapse="-")

#inclusion of 'TSUid' field
sl <- sl[order(sl$PSUid,sl$SSUid),]                                    #Need to be ordered like the result of TSUid below 
PSSUid <- apply(sl[,c("PSUid","SSUid")],1,paste,collapse=".")
PSSUid <- factor(PSSUid,levels=unique(PSSUid))                         #in order to be sure that 'tapply' at the next line will keep the data in order  
tsuidFun <- function(x) as.numeric(as.character(factor(x,levels=unique(x),labels=1:length(unique(x)))))  #function to be applied to catchCat*landCat*commCatScl*commCat crossed values 
                                                                                                         #by PSSUid values to create TSUid field in sl          
ccat <- apply(sl[,c("catchCat","landCat","commCatScl","commCat")],1,paste,collapse="")
sl$TSUid <- unlist(tapply(ccat,list(PSSUid),tsuidFun))



#-------------------------------------------------------------------------------
# Addition of the fields PSUid, SSUid, TSUid, sort, time, space and technical in hl 
# from sl
#-------------------------------------------------------------------------------

hl <- merge(HL,sl[,c("sampType",
                     "landCtry",
                     "vslFlgCtry",
                     "year",
                     "proj",
                     "trpCode",
                     "staNum",
                     "spp",
                     "catchCat",
                     "landCat",
                     "commCatScl",
                     "commCat",
                     "subSampCat",
                     "sort",
                     "time",
                     "space",
                     "technical",
                     "PSUid",
                     "SSUid",
                     "TSUid")],sort=FALSE,all.x=TRUE)

hl <- hl[order(hl$PSUid,hl$SSUid,hl$TSUid),]




#-------------------------------------------------------------------------------
# Addition of the fields PSUid, SSUid, time, space and technical in ca 
# from hh
#-------------------------------------------------------------------------------
#if tech strata = commercial cat, technical field must be included at first
if (cc) CA$technical <- CA$commCat
#Addition of the fields PSUid, SSUid, TSUid,time, space and technical in ca from hh table --> presence of NAs will mean that information is only linked to tr table 
ca <- merge(CA,HH[,c("sampType",
                     "landCtry",
                     "vslFlgCtry",
                     "year",
                     "proj",
                     "trpCode",
                     "staNum",
                     "time",
                     "space",
                     "technical",
                     "PSUid",
                     "SSUid")],sort=FALSE,all.x=TRUE)

# For ca table, there are 2 cases to distinguish : the part of ca that is linked to hh table,
# and the part of ca that is only linked to tr

#the ca part unlinked to hh is isolated as caU...
caU <- ca[is.na(ca$PSUid),]
#...and we keep the other part as ca
ca <- ca[!is.na(ca$PSUid),]


    #-------------------------------------------------------------------------------
    # Creation of the 3 stratification fields in caU (time, space and technical)
    #-------------------------------------------------------------------------------

#'semester' field is added to caU table
caU$semester <- ceiling(caU$quarter/2)

if (is.na(timeStrata)) {
  caU$time <- NA 
} else {
  caU$time <- caU[,timeStrata]}
  
if (is.na(spaceStrata)){ 
  caU$space <- NA 
} else {
  caU$space <- caU[,spaceStrata]}

#technical case is a bit different because only "commCat" values are in ca table  
if (cc){ 
  caU$technical <- caU$commCat
} else {
  caU$technical <- NA } 
    
    
    #-------------------------------------------------------------------------------
    # Addition of PSUid and sort fields in caU table
    #-------------------------------------------------------------------------------

# PSUid related to ca unlinked to hh are currently NAs in tr
# the new PSUid values must begin where current PSUid values in tr table end 
begin <- max(tr$PSUid,na.rm=TRUE)
#PSUid is defined by the combination of a trip, time, space and technical strata
psuid <- apply(caU[,c("sampType","landCtry","vslFlgCtry","proj","trpCode","time","space","technical")],1,paste,collapse=":-:")
caU$PSUid <- as.numeric(as.character(factor(psuid,levels=unique(psuid),labels=(1:length(unique(psuid)))+begin)))
#we can now paste 'ca' and 'caU'... 
neoca <- rbind.data.frame(ca,caU[,names(ca)])           # <--- caU gained the field "semester" above that needs to be removed
#...and creating 'sort' field
fields <- c("catchCat","landCat","commCatScl","commCat")
neoca$sort <- apply(neoca[,fields],1,paste,collapse="-")



    #-------------------------------------------------------------------------------
    # Inclusion of the fields PSUid, time, space and technical in tr table from caU table
    #-------------------------------------------------------------------------------

colN <- c("sampType","landCtry","vslFlgCtry","year","proj","trpCode")    
index <- match(apply(tr[,colN],1,paste,collapse=""),apply(caU[,colN],1,paste,collapse=""))
indexTr <- (1:nrow(tr))[!is.na(index)]
#index & indexTr are linking tr and caU datas
tr$PSUid[indexTr] <- caU$PSUid[index[!is.na(index)]]
tr$time[indexTr] <- caU$time[index[!is.na(index)]]
tr$space[indexTr] <- caU$space[index[!is.na(index)]]
#no technical strata for tr table part that is only linked to ca
#ordering by PSUid field
tr <- tr[order(tr$PSUid),]


#-------------------------------------------------------------------------------
# Finally, creation of 'csDataCons' object and use of 'coerceCons' function to 
# convert columns following default 'csDataCons' object
#-------------------------------------------------------------------------------

csc <- new("csDataCons")

	tr <- tr[,match(names(tr(csc)),names(tr))]
  rownames(tr) <- 1:nrow(tr)
	
  hh <- HH[,match(names(hh(csc)),names(HH))]
  rownames(hh) <- 1:nrow(hh)
  
	sl <- sl[,match(names(sl(csc)),names(sl))]
  rownames(sl) <- 1:nrow(sl)
  
	hl <- hl[,match(names(hl(csc)),names(hl))]
  rownames(hl) <- 1:nrow(hl)
  
	ca <- neoca[,match(names(ca(csc)),names(neoca))]
  rownames(ca) <- 1:nrow(ca)  
  
new("csDataCons",desc=desc,tr=coerceCons(tr,csc@tr),hh=coerceCons(hh,csc@hh),sl=coerceCons(sl,csc@sl),hl=coerceCons(hl,csc@hl),ca=coerceCons(ca,csc@ca))
})





setMethod("csDataCons", signature("csDataVal","missing"), function(object,desc="Unknown stock", ...){

	csDataCons(object,strIni(),desc=desc,...)
})

setMethod("csDataCons", signature("missing","missing"), function(desc="Unknown stock", ...){

	new("csDataCons", desc=desc)
})	





#====================================================================
# Accessor functions
#====================================================================

setMethod("tr", signature("csDataCons"), function(object, ...){
	object@tr
	}
)

setMethod("hh", signature("csDataCons"), function(object, ...){
	object@hh
	}
)

setMethod("sl", signature("csDataCons"), function(object, ...){
	object@sl
	}
)

setMethod("hl", signature("csDataCons"), function(object, ...){
	object@hl
	}
)

setMethod("ca", signature("csDataCons"), function(object, ...){
	object@ca
	}
)

setMethod("desc", signature("csDataCons"), function(object, ...){
	object@desc
	}
)

#====================================================================
# 'Head' and 'Tail' functions
#====================================================================

setMethod("head", signature("csDataCons"), function(x, ...){
  object <- new("csDataCons",desc=x@desc)
  object@tr <- head(x@tr)
  object@hh <- head(x@hh)
  object@sl <- head(x@sl)
  object@hl <- head(x@hl)
  object@ca <- head(x@ca)
  return(object)  
	}
)

setMethod("tail", signature("csDataCons"), function(x, ...){
  object <- new("csDataCons",desc=x@desc)
  object@tr <- tail(x@tr)
  object@hh <- tail(x@hh)
  object@sl <- tail(x@sl)
  object@hl <- tail(x@hl)
  object@ca <- tail(x@ca)
  return(object)  
	}
)

#====================================================================
# 'summary' function
#====================================================================

setMethod("summary", signature("csDataCons"), function(object, ...){
  ll <- list()
  ll$desc <- object@desc
  ll$tr <- summary(object@tr)
  ll$hh <- summary(object@hh)
  ll$sl <- summary(object@sl)
  ll$hl <- summary(object@hl)
  ll$ca <- summary(object@ca)
  return(ll)  
	}
)

#====================================================================
# 'dim' function
#====================================================================

setMethod("dim", signature("csDataCons"), function(x){
  ll <- list()
  ll$tr <- dim(x@tr)
  ll$hh <- dim(x@hh)
  ll$sl <- dim(x@sl)
  ll$hl <- dim(x@hl)
  ll$ca <- dim(x@ca)
  return(ll)  
	}
)

#====================================================================
# 'is.' function
#====================================================================

setGeneric("is.csDataCons", function(object){
	standardGeneric("is.csDataCons")
})


setMethod("is.csDataCons","ANY", function(object){
	return(is(object, "csDataCons"))
})

#====================================================================
# rbind2
#====================================================================

setMethod("rbind2", signature(x="csDataCons", y="csDataCons"), function(x,y){

	# get info
	tr1 <- tr(x)
	hh1 <- hh(x)
	sl1 <- sl(x)
	hl1 <- hl(x)
	ca1 <- ca(x)

	tr2 <- tr(y)
	hh2 <- hh(y)
	sl2 <- sl(y)
	hl2 <- hl(y)
	ca2 <- ca(y)

	# bind
	tr <- rbind2(tr1,tr2)
	hh <- rbind2(hh1,hh2)
	sl <- rbind2(sl1,sl2)
	hl <- rbind2(hl1,hl2)
	ca <- rbind2(ca1,ca2)

	# new object
	new("csDataCons", tr=unique(tr), hh=unique(hh), sl=unique(sl), hl=unique(hl), ca=unique(ca))
})

#====================================================================
# subset
#====================================================================

setMethod("subset", signature(x="csDataCons"), function(x,subset,..., table="tr"){

	# get idx
	trpk <- tr(x)$PSUid
	hhfk <- hh(x)$PSUid
	hhpk <- hh(x)$SSUid
	slfk <- sl(x)$SSUid
	slpk <- sl(x)$TSUid
	hlfk <- hl(x)$TSUid
	cafk <- ca(x)$PSUid
	cafk2 <- ca(x)$SSUid
	
	# new idx
	e <- substitute(subset)
	df0 <- do.call(table, list(object=x))
	r <- eval(e, df0, parent.frame())
	
	# subset
	if(table=="tr"){
		tr <- df0[r,]
		hh <- hh[hh$PSUid %in% tr$PSUid]
		sl <- sl[sl$SSUid %in% hh$SSUid]
		hl <- hl[hl$TSUid %in% sl$TSUid]
		ca <- ca[ca$PSUid %in% tr$PSUid]
	} else if (table=="hh"){
		hh <- df0[r,]
		tr <- tr[tr$PSUid %in% unique(hh$PSUid)]
		sl <- sl[sl$SSUid %in% hh$SSUid]
		hl <- hl[hl$TSUid %in% sl$TSUid]
		ca <- ca[ca$PSUid %in% tr$PSUid]
	} else if(table=="sl"){
		sl <- df0[r,]
		tr <- tr[tr$PSUid %in% unique(sl$PSUid)]
		hh <- hh[hh$SSUid %in% unique(sl$SSUid)]
		hl <- hl[hl$TSUid %in% sl$TSUid]
		ca <- ca[ca$PSUid %in% tr$PSUid]
	} else if(table=="hl"){
		hl <- df0[r,]
		tr <- tr[tr$PSUid %in% unique(hl$PSUid)]
		hh <- hh[hh$SSUid %in% unique(hl$SSUid)]
		sl <- sl[sl$TSUid %in% unique(hl$TSUid)]
		ca <- ca[ca$PSUid %in% tr$PSUid]
	} else if(table=="ca"){
		ca <- df0[r,]
		tr <- tr[tr$PSUid %in% unique(ca$PSUid)]
		hh <- hh[hh$PSUid %in% tr$PSUid]
		sl <- sl[sl$SSUid %in% hh$SSUid]
		hl <- hl[hl$TSUid %in% sl$TSUid]
	}		

	# output
	if(nrow(tr)<1) csData()
	if(nrow(hh)<1) csData(tr=tr)
	if(nrow(sl)<1) csData(tr=tr, hh=hh)
	if(nrow(hl)<1) csData(tr=tr, hh=hh, sl=sl)
	if(nrow(ca)<1) csData(tr=tr, hh=hh, sl=sl, hl=hl)
	else csData(tr=tr, hh=hh, sl=sl, hl=hl, ca=ca)
})


#====================================================================
# replacement
#====================================================================
#
#setMethod("csDataCons", signature("csDataVal"), function(object, desc="Unknown stock", ...){
#
#	#------------------------------------------------------------------------------
#	# create SUid
#	#------------------------------------------------------------------------------
#	suid <- createSUid(object)
#
#	#------------------------------------------------------------------------------
#	# time
#	#------------------------------------------------------------------------------
#	# hh
#	hh <- suid$hh
#	tm <- as.Date(hh$date)
#	tm <- paste(hh$year,quarters(tm),sep=".")
#	hh$time <- tm
#	# tr
#	tm <- unique(suid$hh[,c("PSUid", "time")])
#	tr <- merge(suid$tr, tm, all.x=TRUE)
#	# sl
#	sl <- suid$sl
#	sl$time <- hh$time[match(sl$SSUid,hh$SSUid)]
#	# hl
#	hl <- suid$hl
#	hl$time <- hh$time[match(hl$SSUid,hh$SSUid)]
#	# ca
#	ca <- suid$ca
#	ca$time <- paste(ca$year, paste("Q", ca$quarter, sep=""), sep=".")
#	
#	#------------------------------------------------------------------------------
#	# tech
#	#------------------------------------------------------------------------------
#	te1 <- apply(tr[,c("vslLen", "vslPwr", "vslSize", "vslType")], 1,paste, collapse=".")
#	te2 <- apply(hh[,c("foCatNat","foCatEu5","foCatEu6","gear","meshSize","selDev","meshSizeSelDev")], 1,paste, collapse=".")
#	#hh
#	hh$technical <- paste(te1[match(hh$PSUid,tr$PSUid)], te2, sep="-")
#	# tr
#	te <- unique(hh[,c("PSUid", "technical")])
#	tr <- merge(tr, te, all.x=TRUE)
#	# sl
#	sl$technical <- hh$technical[match(sl$SSUid,hh$SSUid)]
#	# hl
#	hl$technical <- hh$technical[match(hl$SSUid,hh$SSUid)]
#	# ca
#	ca$technical <- hh$technical[match(ca$SSUid,hh$SSUid)]
#	
#	# Not relevant
#	
#	#------------------------------------------------------------------------------
#	# space
#	#------------------------------------------------------------------------------
#	# hh
#	hh$space <- apply(hh[,c("area","rect","waterDep")], 1,paste, collapse=".") 
#	# tr
#	# not relevant
#	# sl
#	sl$space <- hh$space[match(sl$SSUid,hh$SSUid)]
#	# hl
#	hl$space <- hh$space[match(hl$SSUid,hh$SSUid)]
#	# ca
#	ca$space <- apply(ca[,c("area","rect")], 1,paste, collapse=".")
#	
#	#------------------------------------------------------------------------------
#	# sort
#	#------------------------------------------------------------------------------
#	# sl
#	sl$sort <- apply(sl[,c("catchCat","landCat","commCatScl", "commCat", "subSampCat")], 1,paste, collapse=".")
#	# hl
#	hl$sort <- sl$sort[match(hl$TSUid,sl$TSUid)]
#	# ca
#	ca$sort <- apply(ca[,c("catchCat","landCat","commCatScl", "commCat")], 1,paste, collapse=".")
#	
#	#------------------------------------------------------------------------------
#	# create csDataCons
#	#------------------------------------------------------------------------------
#	csc <- csDataCons()
#	tr <- tr[,match(names(tr(csc)),names(tr))]
#	hh <- hh[,match(names(hh(csc)),names(hh))]
#	sl <- sl[,match(names(sl(csc)),names(sl))]
#	hl <- hl[,match(names(hl(csc)),names(hl))]
#	ca <- ca[,match(names(ca(csc)),names(ca))]
#	new("csDataCons", tr=tr, hh=hh, sl=sl, hl=hl, ca=ca)
#})
#
#