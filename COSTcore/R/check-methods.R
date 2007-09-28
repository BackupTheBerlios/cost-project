#====================================================================
#
# EJ, 24/09/2007
# check methods to help on validity and constructors
#
#====================================================================

#====================================================================
# A set of methods to check names in tables 
#====================================================================

setGeneric("checkNms", function(object, nms, ...){
	standardGeneric("checkNms")
	}
)

setMethod("checkNms", signature("data.frame", "character"), function(object, nms, ...){
	lnms <- length(nms)
	nms1 <- names(object)
	if(length(nms1)!=lnms) stop(return(FALSE))
#	if(sum(nms1 %in% nms)!=lnms) stop(return(FALSE))
	if(sum(nms1==nms)!=lnms) stop(return(FALSE))
	return(TRUE)
})

# TR
setGeneric("checkTRnms", function(object, ...){
	standardGeneric("checkTRnms")
	}
)

setMethod("checkTRnms", signature(object="data.frame"), function(object, ...){
	rnms <- c("RECORD_TYPE", "SAMPLING_TYPE", "LANDING_COUNTRY", "VESSEL_FLAG_COUNTRY", "YEAR", "PROJECT", "TRIP_NUMBER", "VESSEL_LENGTH", "VESSEL_POWER", "VESSEL_SIZE", "VESSEL_TYPE", "NUMBER_HAULS_SETS", "DAYS_AT_SEA", "VESSEL_ID", "SAMPLING_COUNTRY", "SAMPLING_METHOD")
	if(checkNms(object, rnms)==FALSE) stop("Check table \"tr\" columns' size and names.")
	return(TRUE)
})

# HH
setGeneric("checkHHnms", function(object, ...){
	standardGeneric("checkHHnms")
	}
)

setMethod("checkHHnms", signature(object="data.frame"), function(object, ...){
	nms <- names(object)
	rnms <- c("RECORD_TYPE", "SAMPLING_TYPE", "LANDING_COUNTRY", "VESSEL_FLAG_COUNTRY", "YEAR", "PROJECT", "TRIP_NUMBER", "STATION_NUMBER", "HAUL_VALIDITY", "AGGREGATION_LEVEL", "DATE", "TIME_SHOT", "FISHING_TIME", "POS_START_LAT_DEC", "POS_START_LON_DEC", "POS_STOP_LAT_DEC", "POS_STOP_LON_DEC", "AREA", "RECTANGLE", "MAIN_FISHING_DEPTH", "MAIN_WATER_DEPTH", "FISHING_ACTIVITY_NAT", "FISHING_ACTIVITY_EU_L5", "FISHING_ACTIVITY_EU_L6", "GEAR", "MESH_SIZE", "SELECTION_DEVICE", "MESH_SIZE_IN_SEL_DEV")
	if(checkNms(object, rnms)==FALSE) stop("Check table \"hh\" columns' size and names.")
	return(TRUE)
})

# SL
setGeneric("checkSLnms", function(object, ...){
	standardGeneric("checkSLnms")
	}
)

setMethod("checkSLnms", signature(object="data.frame"), function(object, ...){
	nms <- names(object)
	rnms <- c("RECORD_TYPE", "SAMPLING_TYPE", "LANDING_COUNTRY", "VESSEL_FLAG_COUNTRY", "YEAR", "PROJECT", "TRIP_NUMBER", "STATION_NUMBER", "SPECIES_CODE", "SEX", "CATCH_CATEGORY", "LANDING_CATEGORY", "COMM_SIZE_CAT_SCALE", "COMM_SIZE_CAT", "SUBSAMPLING_CATEGORY", "VALIDITY_CODE", "WEIGHT", "SUBSAMPLE_WEIGHT", "LENGTH_CODE")
	if(checkNms(object, rnms)==FALSE) stop("Check table \"sl\" columns' size and names.")
	return(TRUE)
})

# HL
setGeneric("checkHLnms", function(object, ...){
	standardGeneric("checkHLnms")
	}
)

setMethod("checkHLnms", signature(object="data.frame"), function(object, ...){
	nms <- names(object)
	rnms <- c("RECORD_TYPE", "SAMPLING_TYPE", "LANDING_COUNTRY", "VESSEL_FLAG_COUNTRY", "YEAR", "PROJECT", "TRIP_NUMBER", "STATION_NUMBER", "SPECIES_CODE", "SEX", "CATCH_CATEGORY", "LANDING_CATEGORY", "COMM_SIZE_CAT_SCALE", "COMM_SIZE_CAT", "SUBSAMPLING_CATEGORY", "LENGTH_CLASS", "NUMBER_AT_LENGTH")
	if(checkNms(object, rnms)==FALSE) stop("Check table \"hl\" columns' size and names.")
	return(TRUE)
})

# CA

#====================================================================
# A set of methods to check PK in tables 
#====================================================================

# TR
setGeneric("checkTRpk", function(object, ...){
	standardGeneric("checkTRpk")
	}
)

setMethod("checkTRpk", signature(object="data.frame"), function(object, ...){
	all.equal(object[,1:6], unique(object[,1:6]))
})

# HH
setGeneric("checkHHpk", function(object, ...){
	standardGeneric("checkHHpk")
	}
)

setMethod("checkHHpk", signature(object="data.frame"), function(object, ...){
	all.equal(object[,1:7], unique(object[,1:7]))
})

# SL
setGeneric("checkSLpk", function(object, ...){
	standardGeneric("checkSLpk")
	}
)

setMethod("checkSLpk", signature(object="data.frame"), function(object, ...){
	all.equal(object[,1:14], unique(object[,1:14]))
})

# HL
setGeneric("checkHLpk", function(object, ...){
	standardGeneric("checkHLpk")
	}
)

setMethod("checkHLpk", signature(object="data.frame"), function(object, ...){
	all.equal(object[,1:15], unique(object[,1:15]))
})

# CA
setGeneric("checkCApk", function(object, ...){
	standardGeneric("checkCApk")
	}
)

setMethod("checkCApk", signature(object="data.frame"), function(object, ...){
	all.equal(object[,c(1:16,18,19,27)], unique(object[,c(1:16,18,19,27)]))
})

#====================================================================
# A set of methods to check data integrity 
#====================================================================

# check data integrity
setGeneric("checkDataIntegrity", function(target, current, ...){
	standardGeneric("checkDataIntegrity")
	}
)

setMethod("checkDataIntegrity", signature(target="data.frame", current="data.frame"), function(target, current, ...){
	trg <- apply(target, 2, paste, collapse="")
	crr <- apply(unique(current), 2, paste, collapse="")
	sum(crr %in% trg)==length(crr)
})


# TR
setGeneric("checkHH2TR", function(hh, tr, ...){
	standardGeneric("checkHH2TR")
	}
)

setMethod("checkHH2TR", signature(tr="data.frame", hh="data.frame"), function(tr, hh, ...){
	# coerce to character to avoid problems with factor definition
	tr <- apply(tr[,1:6], 2, as.character)
	hh <- apply(unique(hh[,1:6]),2,as.character)
	all.equal(tr, hh)
})

# HH
setGeneric("checkHH2SL", function(hh, sl, ...){
	standardGeneric("checkHH2SL")
	}
)

setMethod("checkHH2SL", signature(hh="data.frame", sl="data.frame"), function(hh, sl, ...){
	hh <- apply(hh[,1:7], 2, as.character)
	sl <- apply(unique(sl[,1:7]), 2, as.character)
	all.equal(hh, sl)
})

# SL
setGeneric("checkSLpk", function(object, ...){
	standardGeneric("checkSLpk")
	}
)

setMethod("checkSLpk", signature(object="data.frame"), function(object, ...){
	all.equal(object[,1:14], unique(object[,1:14]))
})

# HL
setGeneric("checkHLpk", function(object, ...){
	standardGeneric("checkHLpk")
	}
)

setMethod("checkHLpk", signature(object="data.frame"), function(object, ...){
	all.equal(object[,1:15], unique(object[,1:15]))
})

# CA
setGeneric("checkCApk", function(object, ...){
	standardGeneric("checkCApk")
	}
)

setMethod("checkCApk", signature(object="data.frame"), function(object, ...){
	all.equal(object[,c(1:16,18,19,27)], unique(object[,c(1:16,18,19,27)]))
})
