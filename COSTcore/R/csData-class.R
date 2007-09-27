#===================================
#
# EJ, 24/09/2007
# csData-class
#
#===================================

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

#===================================
# A set of methods to check names in tables 
#===================================

# TR
setGeneric("checkTRnms", function(object, ...){
	standardGeneric("checkTRnms")
	}
)

setMethod("checkTRnms", signature(object="data.frame"), function(object, ...){
	nms <- names(object)
	rnms <- c("RECORD_TYPE", "SAMPLING_TYPE", "LANDING_COUNTRY", "VESSEL_FLAG_COUNTRY", "YEAR", "PROJECT", "TRIP_NUMBER", "VESSEL_LENGTH", "VESSEL_POWER", "VESSEL_SIZE", "VESSEL_TYPE", "NUMBER_HAULS_SETS", "DAYS_AT_SEA", "VESSEL_ID", "SAMPLING_COUNTRY", "SAMPLING_METHOD")
	lrnms <- length(rnms)
	if(length(nms)!=lrnms) stop("Table \"tr\" must have the following names:\n", paste(rnms, collapse=", "))
	if(sum(nms %in% rnms)!=lrnms) stop("Table \"tr\" must have the following names:\n", paste(rnms, collapse=", "))
})

# HH
setGeneric("checkHHnms", function(object, ...){
	standardGeneric("checkHHnms")
	}
)

setMethod("checkHHnms", signature(object="data.frame"), function(object, ...){
	nms <- names(object)
	rnms <- c("RECORD_TYPE", "SAMPLING_TYPE", "LANDING_COUNTRY", "VESSEL_FLAG_COUNTRY", "YEAR", "PROJECT", "TRIP_NUMBER", "STATION_NUMBER", "HAUL_VALIDITY", "AGGREGATION_LEVEL", "DATE", "TIME_SHOT", "FISHING_TIME", "POS_START_LAT_DEC", "POS_START_LON_DEC", "POS_STOP_LAT_DEC", "POS_STOP_LON_DEC", "AREA", "RECTANGLE", "MAIN_FISHING_DEPTH", "MAIN_WATER_DEPTH", "FISHING_ACTIVITY_NAT", "FISHING_ACTIVITY_EU_L5", "FISHING_ACTIVITY_EU_L6", "GEAR", "MESH_SIZE", "SELECTION_DEVICE", "MESH_SIZE_IN_SEL_DEV")
	lrnms <- length(rnms)
	if(length(nms)!=lrnms) stop("Table \"hh\" must have the following names:\n", paste(rnms, collapse=", "))
	if(sum(nms %in% rnms)!=lrnms) stop("Table \"hh\" must have the following names:\n", paste(rnms, collapse=", "))
})

# SL
setGeneric("checkSLnms", function(object, ...){
	standardGeneric("checkSLnms")
	}
)

setMethod("checkSLnms", signature(object="data.frame"), function(object, ...){
	nms <- names(object)
	rnms <- c("RECORD_TYPE", "SAMPLING_TYPE", "LANDING_COUNTRY", "VESSEL_FLAG_COUNTRY", "YEAR", "PROJECT", "TRIP_NUMBER", "STATION_NUMBER", "SPECIES_CODE", "SEX", "CATCH_CATEGORY", "LANDING_CATEGORY", "COMM_SIZE_CAT_SCALE", "COMM_SIZE_CAT", "SUBSAMPLING_CATEGORY", "VALIDITY_CODE", "WEIGHT", "SUBSAMPLE_WEIGHT", "LENGTH_CODE")
	lrnms <- length(rnms)
	if(length(nms)!=lrnms) stop("Table \"sl\" must have the following names:\n", paste(rnms, collapse=", "))
	if(sum(nms %in% rnms)!=lrnms) stop("Table \"sl\" must have the following names:\n", paste(rnms, collapse=", "))
})

# HL
setGeneric("checkHLnms", function(object, ...){
	standardGeneric("checkHLnms")
	}
)

setMethod("checkHLnms", signature(object="data.frame"), function(object, ...){
	nms <- names(object)
	rnms <- c("RECORD_TYPE", "SAMPLING_TYPE", "LANDING_COUNTRY", "VESSEL_FLAG_COUNTRY", "YEAR", "PROJECT", "TRIP_NUMBER", "STATION_NUMBER", "SPECIES_CODE", "SEX", "CATCH_CATEGORY", "LANDING_CATEGORY", "COMM_SIZE_CAT_SCALE", "COMM_SIZE_CAT", "SUBSAMPLING_CATEGORY", "LENGTH_CLASS", "NUMBER_AT_LENGTH")
	lrnms <- length(rnms)
	if(length(nms)!=lrnms) stop("Table \"hl\" must have the following names:\n", paste(rnms, collapse=", "))
	if(sum(nms %in% rnms)!=lrnms) stop("Table \"hl\" must have the following names:\n", paste(rnms, collapse=", "))
})

# constructor
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
