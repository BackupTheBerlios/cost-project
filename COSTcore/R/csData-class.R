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

	# check column types
	tys0 <- lapply(tr0,class)
	if(checkTys(tr, tys0)==FALSE) stop("Column types not correct in slot candidate \"tr\".")
	tys0 <- lapply(hh0,class)
	if(checkTys(hh, tys0)==FALSE) stop("Column types not correct in slot candidate \"hh\".")
	tys0 <- lapply(sl0,class)
	if(checkTys(sl, tys0)==FALSE) stop("Column types not correct in slot candidate \"sl\".")
	tys0 <- lapply(hl0,class)
	if(checkTys(hl, tys0)==FALSE) stop("Column types not correct in slot candidate \"hl\".")
	tys0 <- lapply(ca0,class)
	if(checkTys(ca, tys0)==FALSE) stop("Column types not correct in slot candidate \"ca\".")

	# check data integrity
	if(checkDataIntegrity(tr[,1:6], hh[,1:6])==FALSE) stop("Data integrity problem in table \"hh\". Missing related records in \"tr\".")
	if(checkDataIntegrity(hh[,1:7], sl[,1:7])==FALSE) stop("Data integrity problem in table \"sl\". Missing related records in \"hh\".")
	if(checkDataIntegrity(sl[,1:13], hl[,c(1:8,10:14)])==FALSE) stop("Data integrity problem in table \"hl\". Missing related records in \"sl\".")
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
			sampType=as.character(NA), # PK
			landCtry=as.character(NA), # PK
			vslFlgCtry=as.character(NA), # PK
			year=as.numeric(NA), # PK
			proj=as.character(NA), # PK
			trpCode=as.character(NA), # PK
			vslLen=as.numeric(NA), 
			vslPwr=as.numeric(NA), 
			vslSize=as.numeric(NA), 
			vsType=as.character(NA), 
			foNum=as.numeric(NA), 
			daysAtSea=as.numeric(NA), 
			vslId=as.numeric(NA), 
			sampCtry=as.character(NA), 
			sampMeth=as.character(NA),
			stringsAsFactors=F),
		hh=data.frame(
			sampType=as.character(NA), # FK
			landCtry=as.character(NA), # FK
			vslFlgCtry=as.character(NA), # FK
			year=as.numeric(NA), # FK
			proj=as.character(NA), # FK
			trpCode=as.character(NA), # FK
			staNum=as.numeric(NA), # PK
			foVal=as.character(NA),
			aggLev=as.character(NA),
			catReg=as.character(NA),
			sppReg=as.character(NA),
			date=as.character(NA),
			time=as.character(NA),
			foDur=as.numeric(NA),
			latIni=as.numeric(NA),
			lonIni=as.numeric(NA),
			latFin=as.numeric(NA),
			lonFin=as.numeric(NA),
			area=as.character(NA),
			rect=as.character(NA),
			foDep=as.numeric(NA),
			waterDep=as.numeric(NA),
			foCatNat=as.character(NA),
			foCatEu5=as.character(NA),
			foCatEu6=as.character(NA),
			gear=as.character(NA),
			meshSize=as.numeric(NA),
			selDev=as.character(NA),
			meshSizeSelDev=as.numeric(NA),
			stringsAsFactors=F),
		sl=data.frame(
			sampType=as.character(NA), # FK
			landCtry=as.character(NA), # FK
			vslFlgCtry=as.character(NA), # FK
			year=as.numeric(NA), # FK
			proj=as.character(NA), # FK
			trpCode=as.character(NA), # FK
			staNum=as.numeric(NA), # FK
			spp=as.character(NA), # PK 
			catchCat=as.character(NA), # PK 
			landCat=as.character(NA), # PK 
			commCatScl=as.character(NA), # PK
			commCat=as.character(NA), # PK
			subSampCat=as.character(NA), # PK
#			valCode=as.character(NA), 
			wt=as.numeric(NA), 
			subSampWt=as.numeric(NA), 
			lenCode=as.character(NA),
			stringsAsFactors=F),
		hl=data.frame(
			sampType=as.character(NA), # FK
			landCtry=as.character(NA), # FK
			vslFlgCtry=as.character(NA), # FK
			year=as.numeric(NA), # FK
			proj=as.character(NA), # FK
			trpCode=as.character(NA), # FK
			staNum=as.numeric(NA), # FK
			spp=as.character(NA), # FK 
			sex=as.character(NA), # PK
			catchCat=as.character(NA), # FK 
			landCat=as.character(NA), # FK 
			commCatScl=as.character(NA), # FK
			commCat=as.character(NA), # FK
			subSampCat=as.character(NA), # FK
			lenCls=as.numeric(NA), # PK
			lenNum=as.numeric(NA),
			stringsAsFactors=F),
		ca=data.frame(
			sampType=as.character(NA), # FK
			landCtry=as.character(NA), # FK
			vslFlgCtry=as.character(NA), # FK
			year=as.numeric(NA), # FK
			proj=as.character(NA), # FK
			trpCode=as.character(NA), # FK
			staNum=as.numeric(NA), # PK
			quarter=as.numeric(NA), # PK
			month=as.numeric(NA), # PK
			spp=as.character(NA), # PK 
			sex=as.character(NA), # PK
			catchCat=as.character(NA), # PK 
			landCat=as.character(NA), # PK 
			commCatScl=as.character(NA), # PK
			commCat=as.character(NA), # PK
			stock=as.character(NA), # PK
			area=as.character(NA), # PK
			rect=as.character(NA), # PK
			lenCls=as.numeric(NA), # PK
			age=as.numeric(NA), # PK
			fishId=as.numeric(NA), # PK
			lenCode=as.character(NA),
			plusGrp=as.character(NA),
			otoWt=as.numeric(NA),
			otoSide=as.character(NA),
			indWt=as.numeric(NA),
			matScale=as.character(NA),
			matStage=as.character(NA),
			stringsAsFactors=F)
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

setMethod("csData", signature("data.frame", "missing", "missing", "missing", "missing"), function(tr, desc="Unknown stock", ...){
	# create object and name columns properly 
	obj <- new("csData")
	
	# tr
	tr0 <- obj@tr
	names(tr) <- names(tr0)
	tr <- coerceDataFrameColumns(tr, tr0)

	# object
	new("csData", tr=tr, desc=desc)
})

setMethod("csData", signature("data.frame", "data.frame", "missing", "missing", "missing"), function(tr, hh, desc="Unknown stock", ...){
	# create object and name columns properly 
	obj <- new("csData")

	# tr
	tr0 <- obj@tr
	names(tr) <- names(tr0)
	tr <- coerceDataFrameColumns(tr, tr0)
	# hh
	hh0 <- obj@hh
	names(hh) <- names(hh0)
	hh <- coerceDataFrameColumns(hh, hh0)

	# object
	new("csData", tr=tr, hh=hh, desc=desc)
})

setMethod("csData", signature("data.frame", "data.frame", "data.frame", "missing", "missing"), function(tr, hh, sl, desc="Unknown stock", ...){
	# create object and name columns properly 
	obj <- new("csData")

	# tr
	tr0 <- obj@tr
	names(tr) <- names(tr0)
	tr <- coerceDataFrameColumns(tr, tr0)
	# hh
	hh0 <- obj@hh
	names(hh) <- names(hh0)
	hh <- coerceDataFrameColumns(hh, hh0)
	# sl
	sl0 <- obj@sl
	names(sl) <- names(sl0)
	sl <- coerceDataFrameColumns(sl, sl0)

	# object
	new("csData", tr=tr, hh=hh, sl=sl, desc=desc)
})

setMethod("csData", signature("data.frame", "data.frame", "data.frame", "data.frame", "missing"), function(tr, hh, sl, hl, desc="Unknown stock", ...){
	# create object and name columns properly 
	obj <- new("csData")

	# tr
	tr0 <- obj@tr
	names(tr) <- names(tr0)
	tr <- coerceDataFrameColumns(tr, tr0)
	# hh
	hh0 <- obj@hh
	names(hh) <- names(hh0)
	hh <- coerceDataFrameColumns(hh, hh0)
	# sl
	sl0 <- obj@sl
	names(sl) <- names(sl0)
	sl <- coerceDataFrameColumns(sl, sl0)
	# hl
	hl0 <- obj@hl
	names(hl) <- names(hl0)
	hl <- coerceDataFrameColumns(hl, hl0)

	new("csData", tr=tr, hh=hh, sl=sl, hl=hl, desc=desc)
})

setMethod("csData", signature("data.frame", "data.frame", "data.frame", "data.frame", "data.frame"), function(tr, hh, sl, hl, ca,  desc="Unknown stock",...){
	# create object and name columns properly 
	obj <- new("csData")
	# tr
	tr0 <- obj@tr
	names(tr) <- names(tr0)
	tr <- coerceDataFrameColumns(tr, tr0)
	# hh
	hh0 <- obj@hh
	names(hh) <- names(hh0)
	hh <- coerceDataFrameColumns(hh, hh0)
	# sl
	sl0 <- obj@sl
	names(sl) <- names(sl0)
	sl <- coerceDataFrameColumns(sl, sl0)
	# hl
	hl0 <- obj@hl
	names(hl) <- names(hl0)
	hl <- coerceDataFrameColumns(hl, hl0)
	# ca
	ca0 <- obj@ca
	names(ca) <- names(ca0)
	ca <- coerceDataFrameColumns(ca, ca0)

	# object
	new("csData", tr=tr, hh=hh, sl=sl, hl=hl, ca=ca, desc=desc)
})

setMethod("csData", signature("data.frame", "missing", "missing", "missing", "data.frame"), function(tr, hh, sl, hl, ca,  desc="Unknown stock",...){
	# create object and name columns properly 
	obj <- new("csData")
	# tr
	tr0 <- obj@tr
	names(tr) <- names(tr0)
	tr <- coerceDataFrameColumns(tr, tr0)
	# ca
	ca0 <- obj@ca
	names(ca) <- names(ca0)
	ca <- coerceDataFrameColumns(ca, ca0)

	# object
	new("csData", tr=tr, ca=ca, desc=desc)
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
	tr <- read.csv(tr,...)
	hh <- read.csv(hh,...)
	sl <- read.csv(sl,...)
	hl <- read.csv(hl,...)

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
	csData(tr=tr, hh=hh, sl=sl, hl=hl, desc=desc)
})


setMethod("csData", signature("character", "character", "character", "character", "character"), function(tr, hh, sl, hl, ca, desc="Unknown stock", ...){

	# read CSV files
	# ToDo
	tr <- read.csv(tr,...)
	hh <- read.csv(hh,...)
	sl <- read.csv(sl,...)
	hl <- read.csv(hl,...)
	ca <- read.csv(ca,...)

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
	csData(tr=tr, hh=hh, sl=sl, hl=hl, ca=ca, desc=desc)
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

#====================================================================
# 'Head' and 'Tail' functions
#====================================================================

setMethod("head", signature("csData"), function(x, ...){
  object <- new("csData",desc=x@desc)
  object@tr <- head(x@tr)
  object@hh <- head(x@hh)
  object@sl <- head(x@sl)
  object@hl <- head(x@hl)
  object@ca <- head(x@ca)
  return(object)  
	}
)

setMethod("tail", signature("csData"), function(x, ...){
  object <- new("csData",desc=x@desc)
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

setMethod("summary", signature("csData"), function(object, ...){
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

setMethod("dim", signature("csData"), function(x){
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

# setGeneric("is.csData", function(object){
# 	standardGeneric("is.csData")
# })


# setMethod("is.csData","ANY", function(object){
# 	return(is(object, "csData"))
# })

#====================================================================
# rbind2
#====================================================================

setMethod("rbind2", signature(x="csData", y="csData"), function(x,y){

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
	csData(tr=unique(tr), hh=unique(hh), sl=unique(sl), hl=unique(hl), ca=unique(ca))
})

#====================================================================
# subset
#====================================================================

setMethod("subset", signature(x="csData"), function(x,subset,..., table="tr"){

	if(table!="tr") stop("Subseting implemented only for slot tr.")
	
	# get idx
	trpk <- tr(x)[,1:6]
	trpk <- apply(trpk,1,paste,collapse="")
	trpk <- gsub("[[:space:]]","",trpk)

	hhfk <- hh(x)[,1:6]
	hhfk <- apply(hhfk,1,paste,collapse="")
	hhfk <- gsub("[[:space:]]","",hhfk)

	hhpk <- hh(x)[,1:7]
	hhpk <- apply(hhpk,1,paste,collapse="")
	hhpk <- gsub("[[:space:]]","",hhpk)

	slfk <- sl(x)[,1:7]
	slfk <- apply(slfk,1,paste,collapse="")
	slfk <- gsub("[[:space:]]","",slfk)

	slpk <- sl(x)[,1:13]
	slpk <- apply(slpk,1,paste,collapse="")
	slpk <- gsub("[[:space:]]","",slpk)

	hlfk <- hl(x)[,c(1:8,10:14)]
	hlfk <- apply(hlfk,1,paste,collapse="")
	hlfk <- gsub("[[:space:]]","",hlfk)

#	hlpk <- hl(x)[,1:15]
#	hlpk <- apply(hlpk,1,paste,collapse="")
#	hlpk <- gsub("[[:space:]]","",hlpk)

	cafk <- ca(x)[,c(1:6)]
	cafk <- apply(cafk,1,paste,collapse="")
	cafk <- gsub("[[:space:]]","",cafk)

#	capk <- ca(x)[,c(1:21)]
#	capk <- apply(capk,1,paste,collapse="")
#	capk <- gsub("[[:space:]]","",capk)
	
	# new idx
	e <- substitute(subset)
	df0 <- do.call(table, list(object=x))
	r <- eval(e, df0, parent.frame())
	
	# subset
	tr <- df0[r,]
	tridx <- apply(tr[,1:6],1,paste,collapse="")
	tridx <- gsub("[[:space:]]","",tridx)
	hh <- hh(x)[hhfk %in% tridx,]	
	hhidx <- apply(hh[,1:7],1,paste,collapse="")
	hhidx <- gsub("[[:space:]]","",hhidx)
	sl <- sl(x)[slfk %in% hhidx,]	
	slidx <- apply(sl[,1:13],1,paste,collapse="")
	slidx <- gsub("[[:space:]]","",slidx)
	hl <- hl(x)[hlfk %in% slidx,]	
	ca <- ca(x)[cafk %in% tridx,]

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

if (!isGeneric("replace")) setGeneric("replace")


