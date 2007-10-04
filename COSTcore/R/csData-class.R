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

valcsData <- function(object){

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
	if(checkTRpk(tr)==FALSE) stop("Primary key not unique in slot candidate \"tr\".")
	if(checkHHpk(hh)==FALSE) stop("Primary key not unique in slot candidate \"hh\".")
	if(checkSLpk(sl)==FALSE) stop("Primary key not unique in slot candidate \"sl\".")
	if(checkHLpk(hl)==FALSE) stop("Primary key not unique in slot candidate \"hl\".")
	if(checkCApk(ca)==FALSE) stop("Primary key not unique in slot candidate \"ca\".")

	# check data integrity
	if(checkDataIntegrity(tr[,1:6], hh[,1:6])==FALSE) stop("Data integrity problem in table \"hh\". Missing related records in \"tr\".")
	if(checkDataIntegrity(hh[,1:7], sl[,1:7])==FALSE) stop("Data integrity problem in table \"sl\". Missing related records in \"hh\".")
	if(checkDataIntegrity(sl[,1:14], hl[,1:14])==FALSE) stop("Data integrity problem in table \"hl\". Missing related records in \"sl\".")
#	if(checkDataIntegrity(tr[,1:6], ca[,1:6])==FALSE) stop("Data integrity problem in table \"ca\". Missing related records in \"tr\".")

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
	),
	validity=valcsData
)

#====================================================================
# Class constructor
#====================================================================
setGeneric("csData", function(tr, hh, sl, hl, ca, ...){
	standardGeneric("csData")
	}
)

setMethod("csData", signature("data.frame", "data.frame", "data.frame", "data.frame", "missing"), function(tr, hh, sl, hl, desc="Unknown stock", ...){
	# create object and name columns properly 
	obj <- new("csData")
	names(tr) <- names(obj@tr)
	names(hh) <- names(obj@hh)
	names(sl) <- names(obj@sl)
	names(hl) <- names(obj@hl)
	new("csData", tr=tr, hh=hh, sl=sl, hl=hl, desc=desc)
})

setMethod("csData", signature("data.frame", "data.frame", "data.frame", "data.frame", "data.frame"), function(tr, hh, sl, hl, ca,  desc="Unknown stock",...){
	# create object and name columns properly 
	obj <- new("csData")
	names(tr) <- names(obj@tr)
	names(hh) <- names(obj@hh)
	names(sl) <- names(obj@sl)
	names(hl) <- names(obj@hl)
	names(ca) <- names(obj@ca)
	new("csData", tr=tr, hh=hh, sl=sl, hl=hl, ca=ca, desc=desc)
})

setMethod("csData", signature("missing", "missing", "missing", "missing", "missing"), function(desc="Unknown stock", ...){
	# create object and name columns properly 
	new("csData", desc=desc)
})

#====================================================================
# IO constructor
#====================================================================

setMethod("csData", signature("character", "character", "character", "character", "missing"), function(tr, hh, sl, hl, desc="Unknown stock", ...){

	# read CSV files
	# ToDo
	tr <- read.csv(tr)
	hh <- read.csv(hh)
	sl <- read.csv(sl)
	hl <- read.csv(hl)

	# check names are correct
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
	new("csData", tr=tr, hh=hh, sl=sl, hl=hl, desc=desc)
})


setMethod("csData", signature("character", "character", "character", "character", "character"), function(tr, hh, sl, hl, ca, desc="Unknown stock", ...){

	# read CSV files
	# ToDo
	tr <- read.csv(tr)
	hh <- read.csv(hh)
	sl <- read.csv(sl)
	hl <- read.csv(hl)
	ca <- read.csv(ca)

	# check names are correct
	checkTRnms(tr)
	checkHHnms(hh)
	checkSLnms(sl)
	checkHLnms(hl)
  checkCAnms(ca)

	# remove record type 
	tr <- tr[,-1]
	hh <- hh[,-1]
	sl <- sl[,-1]
	hl <- hl[,-1]
  ca <- ca[,-1]

	# create object and name columns properly 
	obj <- new("csData")
	names(tr) <- names(obj@tr)
	names(hh) <- names(obj@hh)
	names(sl) <- names(obj@sl)
	names(hl) <- names(obj@hl)
	names(ca) <- names(obj@ca)
	new("csData", tr=tr, hh=hh, sl=sl, hl=hl, ca=ca, desc=desc)
})



#====================================================================
# Accessor functions
#====================================================================

setGeneric("tr", function(object, ...){
	standardGeneric("tr")
	}
)

setMethod("tr", signature("csData"), function(object, ...){
	object@tr
	}
)

setGeneric("hh", function(object, ...){
	standardGeneric("hh")
	}
)

setMethod("hh", signature("csData"), function(object, ...){
	object@hh
	}
)

setGeneric("sl", function(object, ...){
	standardGeneric("sl")
	}
)

setMethod("sl", signature("csData"), function(object, ...){
	object@sl
	}
)

setGeneric("hl", function(object, ...){
	standardGeneric("hl")
	}
)

setMethod("hl", signature("csData"), function(object, ...){
	object@hl
	}
)

setGeneric("ca", function(object, ...){
	standardGeneric("ca")
	}
)

setMethod("ca", signature("csData"), function(object, ...){
	object@ca
	}
)

setMethod("desc", signature("csData"), function(object, ...){
	object@desc
	}
)
