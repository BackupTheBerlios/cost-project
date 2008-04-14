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

tempStrata <- objStrat@tempStrata 
spaceStrata <- objStrat@spaceStrata 
techStrata <- objStrat@techStrata 
sorting <- objStrat@sorting
tpRec <- objStrat@tpRec
spRec <- objStrat@spRec
tcRec <- objStrat@tcRec


TR <- object@tr 
HH <- object@hh 
SL <- object@sl
HL <- object@hl 
CA <- object@ca

#if techStrata="commCat" --> indicator cc
cc <- FALSE
if (!is.null(techStrata)) {
  if (techStrata=="commCat") cc <- TRUE}

if (cc) {
  HH$ind <- 1:nrow(HH)
  HH <- merge(HH,unique(SL[,c("sampType","landCtry","vslFlgCtry","year","proj","trpCode","staNum","commCat")]),sort=FALSE,all=TRUE)
#order
  HH <- HH[order(HH$ind),]}


HH$month <- sapply(HH$date,function(x) as.numeric(strsplit(x,"-")[[1]][2]))
HH$quarter <- ceiling(HH$month/3)
HH$semester <- ceiling(HH$month/6)      

if (is.null(tempStrata)) {
  HH$time <- NA 
  tpRec <- NULL} 
else 
  HH$time <- HH[,tempStrata]  
     
if (is.null(spaceStrata)) {
  HH$space <- NA 
  spRec <- NULL} 
else 
  HH$space <- HH[,spaceStrata]
  
if (is.null(techStrata)) {
  HH$technical <- NA 
  tcRec <- NULL} 
else 
  HH$technical <- HH[,techStrata]


#recoding process
if (!is.null(tpRec)) {
  Typ <- class(HH$time) 
  HH$time <- factor(HH$time) 
  Lev <- levels(HH$time)[!levels(HH$time)%in%tpRec$from]
  HH$time <- factor(HH$time,levels=c(Lev,tpRec$from),labels=c(Lev,tpRec$to))
  eval(parse('',text=paste("HH$time <- as.",Typ,"(as.character(HH$time))",sep="")))}
  
if (!is.null(spRec)) {
  Typ <- class(HH$space) 
  HH$space <- factor(HH$space) 
  Lev <- levels(HH$space)[!levels(HH$space)%in%spRec$from]
  HH$space <- factor(HH$space,levels=c(Lev,spRec$from),labels=c(Lev,spRec$to)) 
  eval(parse('',text=paste("HH$space <- as.",Typ,"(as.character(HH$space))",sep="")))}
  
if (!is.null(tcRec)) {
  Typ <- class(HH$technical) 
  HH$technical <- factor(HH$technical)
  Lev <- levels(HH$technical)[!levels(HH$technical)%in%tcRec$from]
  HH$technical <- factor(HH$technical,levels=c(Lev,tcRec$from),labels=c(Lev,tcRec$to))
  eval(parse('',text=paste("HH$technical <- as.",Typ,"(as.character(HH$technical))",sep="")))}
                        

#PSUid is identified by a trip, a time, space and technical strata
psuid <- apply(HH[,c("sampType","landCtry","vslFlgCtry","proj","trpCode","time","space","technical")],1,paste,collapse="::")
HH$PSUid <- factor(psuid,levels=unique(psuid),labels=1:length(unique(psuid)))
#SSUid
HH$SSUid <- as.numeric(unlist(tapply(HH$PSUid,list(HH$PSUid),function(x) 1:length(x))))




#TR table
tr <- merge(TR,unique(HH[,c("sampType","landCtry","vslFlgCtry","year","proj","trpCode","PSUid","time","space","technical")]),sort=FALSE,all=TRUE)
#tr$foNum & tr$daysAtSea calculation
#only information corresponding to aggLev="H" in HH is concerned 
HHh <- HH[HH$aggLev=="H",]
#object to be inserted in tr...
tab <- spdAgreg(list(nOP=HHh$staNum,nDay=HHh$date),list(PSUid=HHh$PSUid),function(x) length(unique(x))) 
tab <- tab[order(as.numeric(as.character(tab$PSUid))),]
#insert & substitute
tr2 <- merge(tr,tab,sort=FALSE,all.x=TRUE)
tr2$foNum[!is.na(tr2$nOP)] <- tr2$nOP[!is.na(tr2$nOP)]
tr2$daysAtSea[!is.na(tr2$nDay)] <- tr2$nDay[!is.na(tr2$nDay)]
tr <- tr2




#SL table
if (cc) SL$technical <- SL$commCat
sl <- merge(SL,HH[,c("sampType","landCtry","vslFlgCtry","year","proj","trpCode","staNum","time","space","technical","PSUid","SSUid")],sort=FALSE,all.x=TRUE)

if (is.null(sorting)) {
  fields <- NULL
} else {
  sorti <- c(2,4,5)
  names(sorti) <- c("catchCat","commCat","subSampcat")
  fields <- c("catchCat","landCat","commCatScl","commCat","subSampcat")[1:sorti[sorting]]}
  
tsuid <- apply(sl[,c("PSUid","SSUid",fields,"proj","trpCode","time","space","technical")],1,paste,collapse="::")
pss <- apply(sl[,c("PSUid","SSUid")],1,paste,collapse="::")
TS <- tapply(tsuid,list(pss),function(x) as.character(factor(x,levels=unique(x),labels=1:length(unique(x))))) 
sl$TSUid <- factor(unlist(TS[unique(pss)]))
#sampled and measured weights are aggregated
if (is.null(sorting)) 
  sl$sort <- NA 
else 
  sl$sort <- apply(sl[,fields],1,paste,collapse="-")
Sl <- spdAgreg(list(wt=sl$wt,subSampWt=sl$subSampWt),list(sampType=sl$sampType,landCtry=sl$landCtry,vslFlgCtry=sl$vslFlgCtry,year=sl$year,proj=sl$proj,trpCode=sl$trpCode,
                                                             staNum=sl$staNum,spp=sl$spp,lenCode=sl$lenCode,time=sl$time,space=sl$space,technical=sl$technical,sort=sl$sort,PSUid=sl$PSUid,
                                                             SSUid=sl$SSUid,TSUid=sl$TSUid),sum,na.rm=TRUE)
#order
SLl <- Sl[order(as.numeric(as.character(Sl$PSUid)),as.numeric(as.character(Sl$SSUid)),as.numeric(as.character(Sl$TSUid))),]




#HL table
if (cc) HL$technical <- HL$commCat  
if (is.null(sorting)) 
  HL$sort <- NA 
else 
  HL$sort <- apply(HL[,fields],1,paste,collapse="-")
hl <- merge(HL,SLl[,c("sampType","landCtry","vslFlgCtry","year","proj","trpCode","staNum","spp","sort","time","space","technical","PSUid","SSUid","TSUid")],sort=FALSE,all.x=TRUE)
Hl <- spdAgreg(list(lenNum=hl$lenNum),list(sampType=hl$sampType,landCtry=hl$landCtry,vslFlgCtry=hl$vslFlgCtry,year=hl$year,proj=hl$proj,trpCode=hl$trpCode,
                                                             staNum=hl$staNum,spp=hl$spp,sex=hl$sex,lenCls=hl$lenCls,time=hl$time,space=hl$space,technical=hl$technical,sort=hl$sort,
                                                             PSUid=hl$PSUid,SSUid=hl$SSUid,TSUid=hl$TSUid),sum,na.rm=TRUE)
#order
HLl <- Hl[order(as.numeric(as.character(Hl$PSUid)),as.numeric(as.character(Hl$SSUid)),as.numeric(as.character(Hl$TSUid))),]




#CA table
  #part linked to HH
if (cc) CA$technical <- CA$commCat
ca <- merge(CA,HH[,c("sampType","landCtry","vslFlgCtry","year","proj","trpCode","staNum","time","space","technical","PSUid","SSUid")],sort=FALSE,all.x=TRUE)
  #part not linked to HH
ca2 <- ca[is.na(ca$PSUid),]
ca <- ca[!is.na(ca$PSUid),]
#strata fields must be created
ca2$semester <- ceiling(ca2$quarter/2)

if (is.null(tempStrata)) 
  ca2$time <- NA 
else 
  ca2$time <- ca2[,tempStrata]
  
if (is.null(spaceStrata)) 
  ca2$space <- NA 
else 
  ca2$space <- ca2[,spaceStrata]
  
if (is.null(techStrata)) 
  ca2$technical <- NA 
else {
  if (!techStrata%in%names(ca2)) 
    ca2$technical <- NA 
  else 
    ca2$technical <- ca2[,techStrata]}
    
deb <- max(as.numeric(as.character(HH$PSUid)),na.rm=TRUE)
psuid2 <- apply(ca2[,c("sampType","landCtry","vslFlgCtry","proj","trpCode","time","space","technical")],1,paste,collapse="::")
ca2$PSUid <- ca2$psuid <- factor(psuid2,levels=unique(psuid2),labels=(1:length(unique(psuid2)))+deb)
ca2$TIME <- ca2$time
ca2$SPACE <- ca2$space
ca2$TECHNICAL <- ca2$technical  

#insert PSUid in tr
neotr <- merge(tr,unique(ca2[,c("sampType","landCtry","vslFlgCtry","year","proj","trpCode","psuid","TIME","SPACE","TECHNICAL")]),sort=FALSE,all.x=TRUE)
neotr$PSUid <- as.character(neotr$PSUid)
neotr$psuid <- as.character(neotr$psuid)
neotr$PSUid[!is.na(neotr$psuid)] <- neotr$psuid[!is.na(neotr$psuid)]
neotr$time[!is.na(neotr$TIME)] <- neotr$TIME[!is.na(neotr$TIME)]
neotr$space[!is.na(neotr$SPACE)] <- neotr$SPACE[!is.na(neotr$SPACE)]
neotr$technical[!is.na(neotr$TECHNICAL)] <- neotr$TECHNICAL[!is.na(neotr$TECHNICAL)]
#order
neotr <- neotr[order(as.numeric(as.character(neotr$PSUid))),]

neoca <- rbind(ca,ca2[,-(c(-1,0)+ncol(ca2))])
if (is.null(sorting)) 
  neoca$sort <- NA 
else 
  neoca$sort<- apply(neoca[,fields],1,paste,collapse="-")

csc <- new("csDataCons")
	tr <- neotr[,match(names(tr(csc)),names(neotr))] ; rownames(tr) <- 1:nrow(tr)
	hh <- HH[,match(names(hh(csc)),names(HH))] ; rownames(hh) <- 1:nrow(hh)
	sl <- SLl[,match(names(sl(csc)),names(SLl))] ; rownames(sl) <- 1:nrow(sl)
	hl <- HLl[,match(names(hl(csc)),names(HLl))] ; rownames(hl) <- 1:nrow(hl)
	ca <- neoca[,match(names(ca(csc)),names(neoca))] ; rownames(ca) <- 1:nrow(ca)  
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