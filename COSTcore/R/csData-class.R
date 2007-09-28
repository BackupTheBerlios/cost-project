#====================================================================
#
# EJ, 24/09/2007
# csData-class
#
#====================================================================

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


#===================================
#
# Data hierarchy
# tr --> hh --> sl --> hl
# tr --> ca
#  
#===================================

validcsData <- function(object){

	tr <- object@tr
	hh <- object@hh
	sl <- object@sl
	hl <- object@hl
	ca <- object@ca

	# I will rely o the prototype to check col names and size. I'm not sure it's a good strategy !
	obj <- new("csData")
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
	
	# check PK
	if(checkTRpk(tr)==FALSE) stop("Data integrity problem in slot candidate \"tr\". PK not unique.")
	if(checkHHpk(hh)==FALSE) stop("Data integrity problem in slot candidate \"hh\". PK not unique.")
	if(checkSLpk(sl)==FALSE) stop("Data integrity problem in slot candidate \"sl\". PK not unique.")
	if(checkHLpk(hl)==FALSE) stop("Data integrity problem in slot candidate \"hl\". PK not unique.")
	if(checkCApk(ca)==FALSE) stop("Data integrity problem in slot candidate \"ca\". PK not unique.")

	# check data integrity
	if(checkTRpk(tr)==FALSE) stop("Data integrity problem in table \"tr\". PK not unique.")
	if(checkHHpk(hh)==FALSE) stop("Data integrity problem in table \"hh\". PK not unique.")
	if(checkSLpk(sl)==FALSE) stop("Data integrity problem in table \"sl\". PK not unique.")
	if(checkHLpk(hl)==FALSE) stop("Data integrity problem in table \"hl\". PK not unique.")
	if(checkCApk(ca)==FALSE) stop("Data integrity problem in table \"ca\". PK not unique.")

	# Everything is fine
	return(TRUE)
}

setClass("csData",
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
			sampType=NA, # PK
			landCtry=NA, # PK
			vslFlgCtry=NA, # PK
			year=NA, # PK
			proj=NA, # PK
			trpNum=NA, # PK
			vslLen=NA, 
			vslPwr=NA, 
			vslSize=NA, 
			vsType=NA, 
			foNum=NA, 
			daysAtSea=NA, 
			vslId=NA, 
			sampCtry=NA, 
			sampMeth=NA),
		hh=data.frame(
			sampType=NA, # FK
			landCtry=NA, # FK
			vslFlgCtry=NA, # FK
			year=NA, # FK
			proj=NA, # FK
			trpNum=NA, # FK
			staNum=NA, # PK
			foVal=NA,
			aggLev=NA,
			date=NA,
			timeShot=NA,
			foDur=NA,
			latIni=NA,
			lonIni=NA,
			latFin=NA,
			lonFin=NA,
			area=NA,
			rect=NA,
			foDep=NA,
			waterDep=NA,
			foCatNat=NA,
			foCatEu5=NA,
			foCatEu6=NA,
			gear=NA,
			meshSize=NA,
			selDev=NA,
			meshSizeSelDev=NA),
		sl=data.frame(
			sampType=NA, # FK
			landCtry=NA, # FK
			vslFlgCtry=NA, # FK
			year=NA, # FK
			proj=NA, # FK
			trpNum=NA, # FK
			staNum=NA, # FK
			spp=NA, # PK 
			sex=NA, # PK
			catchCat=NA, # PK 
			landCat=NA, # PK 
			commCatScl=NA, # PK
			commCat=NA, # PK
			subSampCat=NA, # PK
			valCode=NA, 
			wt=NA, 
			subSampWt=NA, 
			lenCode=NA),
		hl=data.frame(
			sampType=NA, # FK
			landCtry=NA, # FK
			vslFlgCtry=NA, # FK
			year=NA, # FK
			proj=NA, # FK
			trpNum=NA, # FK
			staNum=NA, # FK
			spp=NA, # PK 
			sex=NA, # PK
			catchCat=NA, # PK 
			landCat=NA, # PK 
			commCatScl=NA, # PK
			commCat=NA, # PK
			subSampCat=NA, # PK
			lenCls=NA, 
			lenNum=NA),
		ca=data.frame(
			sampType=NA, # FK
			landCtry=NA, # FK
			vslFlgCtry=NA, # FK
			year=NA, # FK
			proj=NA, # FK
			trpNum=NA, # FK
			staNum=NA, # PK (optional)
			spp=NA, # PK 
			sex=NA, # PK
			catchCat=NA, # PK 
			landCat=NA, # PK 
			commCatScl=NA, # PK
			commCat=NA, # PK
			stock=NA, # PK
			area=NA, # PK
			rect=NA, # PK
			lenCode=NA, # PK
			lenCls=NA, # PK
			age=NA, # PK
			plusGrp=NA, # PK
			otoWt=NA,
			otoSide=NA,
			indWt=NA,
			matScale=NA,
			matStage=NA,
			num=NA,
			fishId=NA) # PK 
	)
)

# constructor (missing ca)
setGeneric("csData", function(object, ...){
	standardGeneric("csData")
	}
)

setMethod("csData", signature(object="list"), function(object, ...){
	if(sum(names(object) %in% c("tr","hh","sl","hl"))!=4)
		stop("The list must have at least the elements \"tr\", \"hh\", \"sl\" and \"hl\".")
	tr <- object$tr
	hh <- object$hh
	sl <- object$sl
	hl <- object$hl
	
	# checks on names (a bad option that needs to be discussed)	
	checkTRnms(tr)
	checkHHnms(hh)
	checkSLnms(sl)
	checkHLnms(hl)

	# remove record type 
	tr <- tr[,-1]
	hh <- hh[,-1]
	sl <- sl[,-1]
	hl <- hl[,-1]
	
	# create object and name columns properly 
	obj <- new("csData")
	names(tr) <- names(obj@tr)
	names(hh) <- names(obj@hh)
	names(sl) <- names(obj@sl)
	names(hl) <- names(obj@hl)
	new("csData", tr=tr, hh=hh, sl=sl, hl=hl)
})

#setMethod("csData", signature("data.frame"), function(tr, hh, sl, hl, ...){
#
#	# checks on names (a bad option that needs to be discussed)	
#	checkTRnms(tr)
#	checkHHnms(hh)
#	checkSLnms(sl)
#	checkHLnms(hl)
#
#	# remove record type 
#	tr <- tr[,-1]
#	hh <- hh[,-1]
#	sl <- sl[,-1]
#	hl <- hl[,-1]
#	
#	# create object and name columns properly 
#	obj <- new("csData")
#	names(tr) <- names(obj@tr)
#	names(hh) <- names(obj@hh)
#	names(sl) <- names(obj@sl)
#	names(hl) <- names(obj@hl)
#	new("csData", tr=tr, hh=hh, sl=sl, hl=hl)
#})
#
