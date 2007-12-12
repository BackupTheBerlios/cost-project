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
			PSUid=NA, # field to be created by the EDA process, one per each record = PK
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
			vslId=NA, # PK
			sampCtry=NA, 
			sampMeth=NA),
		hh=data.frame(
			PSUid=NA, # FK
			staNum=NA, # PK
			SSUid=NA, # field to be created by the EDA process, one per each record = FK+PK
			time=NA,
			space=NA,
			technical=NA,
			aggLev=NA,
			foDur=NA),
		sl=data.frame(
			PSUid=NA, # This field helps on linking with TR table
			SSUid=NA, # FK 
			spp=NA, # PK 
			catchCat=NA, # PK 
			landCat=NA, # PK
			commCatScl=NA, # PK
			commCat=NA, # PK
			subSampCat=NA, # PK
			TSUid=NA, # field to be created by the EDA process, one per each record = FK+PK
			time=NA,
			space=NA,
			technical=NA,
			sort=NA, 
			wt=NA, 
			subSampWt=NA, 
			lenCode=NA),
		hl=data.frame(
			PSUid=NA, # This field helps on linking with TR table
			SSUid=NA, # This field helps on linking with HH table 
			TSUid=NA, # FK
			time=NA,
			space=NA,
			technical=NA,
			sort=NA, 
			lenCls=NA, # PK
			lenNum=NA),
		ca=data.frame(
			PSUid=NA, # FK
			SSUid=NA, # must match SSUid in HH so that info about tech can be used if necessary
			space=NA,
			time=NA,
			technical=NA,
			sort=NA,
			spp=NA, # PK 
			sex=NA, # PK
			catchCat=NA, # PK 
			landCat=NA, # PK 
			commCatScl=NA, # PK
			commCat=NA, # PK
			stock=NA, # PK
			area=NA, # PK
			rect=NA, # PK
			lenCls=NA, # PK
			age=NA, # PK
			fishId=NA, # PK
			lenCode=NA,
			plusGrp=NA,
			otoWt=NA,
			otoSide=NA,
			indWt=NA,
			matScale=NA,
			matStage=NA)
	),
	validity=valcscData
)

#====================================================================
# Class constructor
#====================================================================
setGeneric("csDataCons", function(object, ...){
	standardGeneric("csDataCons")
	}
)

setMethod("csDataCons", signature("csDataVal"), function(object, desc="Unknown stock", ...){
	# create object and name columns properly 
	stop("Not implemented yet !")
})

setMethod("csDataCons", signature("missing"), function(desc="Unknown stock", ...){
	# create object and name columns properly 
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

