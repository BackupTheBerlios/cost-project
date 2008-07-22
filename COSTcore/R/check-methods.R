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
	nms <- toupper(nms)
	lnms <- length(nms)
	nms1 <- toupper(names(object))
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
	rnms <- c("RECORD_TYPE", "SAMPLING_TYPE", "LANDING_COUNTRY", "VESSEL_FLAG_COUNTRY", "YEAR", "PROJECT", "TRIP_NUMBER", "HARBOUR", "VESSEL_LENGTH", "VESSEL_POWER", "VESSEL_SIZE", "VESSEL_TYPE", "NUMBER_HAULS_SETS", "DAYS_AT_SEA", "VESSEL_ID", "SAMPLING_COUNTRY", "SAMPLING_METHOD")
	if(checkNms(object, rnms)==FALSE) stop("Check table \"TR\" columns' size and names.")
	return(TRUE)
})

# HH
setGeneric("checkHHnms", function(object, ...){
	standardGeneric("checkHHnms")
	}
)

setMethod("checkHHnms", signature(object="data.frame"), function(object, ...){
	nms <- names(object)
	rnms <- c("RECORD_TYPE", "SAMPLING_TYPE", "LANDING_COUNTRY", "VESSEL_FLAG_COUNTRY", "YEAR", "PROJECT", "TRIP_NUMBER", "STATION_NUMBER", "FISHING_VALIDITY", "AGGREGATION_LEVEL", "CATCH_REGISTRATION", "SPECIES_REGISTRATION", "DATE", "TIME", "FISHING_DURATION", "POS_START_LAT_DEC", "POS_START_LON_DEC", "POS_STOP_LAT_DEC", "POS_STOP_LON_DEC", "AREA", "RECTANGLE", "SUB_RECTANGLE", "MAIN_FISHING_DEPTH", "MAIN_WATER_DEPTH", "FISHING_ACTIVITY_NAT", "FISHING_ACTIVITY_EU_L5", "FISHING_ACTIVITY_EU_L6", "MESH_SIZE", "SELECTION_DEVICE", "MESH_SIZE_IN_SEL_DEV")
	if(checkNms(object, rnms)==FALSE) stop("Check table \"HH\" columns' size and names.")
	return(TRUE)
})

# SL
setGeneric("checkSLnms", function(object, ...){
	standardGeneric("checkSLnms")
	}
)

setMethod("checkSLnms", signature(object="data.frame"), function(object, ...){
	nms <- names(object)
	rnms <- c("RECORD_TYPE", "SAMPLING_TYPE", "LANDING_COUNTRY", "VESSEL_FLAG_COUNTRY", "YEAR", "PROJECT", "TRIP_NUMBER", "STATION_NUMBER", "SPECIES_CODE", "CATCH_CATEGORY", "LANDING_CATEGORY", "COMM_SIZE_CAT_SCALE", "COMM_SIZE_CAT", "SUBSAMPLING_CATEGORY", "SEX", "WEIGHT", "SUBSAMPLE_WEIGHT", "LENGTH_CODE")
	if(checkNms(object, rnms)==FALSE) stop("Check table \"SL\" columns' size and names.")
	return(TRUE)
})

# HL
setGeneric("checkHLnms", function(object, ...){
	standardGeneric("checkHLnms")
	}
)

setMethod("checkHLnms", signature(object="data.frame"), function(object, ...){
	nms <- names(object)
	rnms <- c("RECORD_TYPE", "SAMPLING_TYPE", "LANDING_COUNTRY", "VESSEL_FLAG_COUNTRY", "YEAR", "PROJECT", "TRIP_NUMBER", "STATION_NUMBER", "SPECIES_CODE", "CATCH_CATEGORY", "LANDING_CATEGORY", "COMM_SIZE_CAT_SCALE", "COMM_SIZE_CAT", "SUBSAMPLING_CATEGORY", "SEX", "LENGTH_CLASS", "NUMBER_AT_LENGTH")
	if(checkNms(object, rnms)==FALSE) stop("Check table \"HL\" columns' size and names.")
	return(TRUE)
})

# CA
setGeneric("checkCAnms", function(object, ...){
	standardGeneric("checkCAnms")
	}
)

setMethod("checkCAnms", signature(object="data.frame"), function(object, ...){
	nms <- names(object)
	rnms <- c("RECORD_TYPE", "SAMPLING_TYPE", "LANDING_COUNTRY", "VESSEL_FLAG_COUNTRY", "YEAR", "PROJECT", "TRIP_NUMBER", "STATION_NUMBER", "QUARTER", "MONTH", "SPECIES_CODE", "SEX", "CATCH_CATEGORY", "LANDING_CATEGORY", "COMM_SIZE_CAT_SCALE", "COMM_SIZE_CAT", "STOCK", "AREA", "RECTANGLE", "SUB_RECTANGLE", "LENGTH_CLASS", "AGE", "SINGLE_FISH_NB", "LENGTH_CODE", "AGE_PLUS_GROUP", "AGEING_METHOD", "OTOLITH_WEIGHT", "OTOLITH_SIDE", "INDIVIDUAL_WEIGHT", "MATURITY_SCALE", "MATURITY_STAGE", "MATURITY_STAGING_METHOD")
	if(checkNms(object, rnms)==FALSE) stop("Check table \"CA\" columns' size and names.")
	return(TRUE)
})

# CL
setGeneric("checkCLnms", function(object, ...){
	standardGeneric("checkCLnms")
	}
)

setMethod("checkCLnms", signature(object="data.frame"), function(object, ...){
	nms <- names(object)
	rnms <- c("RECORD_TYPE", "LANDING_COUNTRY", "VESSEL_FLAG_COUNTRY", "YEAR", "QUARTER", "MONTH", "AREA", "RECTANGLE", "SUB_RECTANGLE", "TAXON", "LANDING_CATEGORY", "COMM_SIZE_CAT_SCALE", "COMM_SIZE_CAT", "FISHING_ACTIVITY_NAT", "FISHING_ACTIVITY_EU_L5", "FISHING_ACTIVITY_EU_L6", "HARBOUR", "UNALL_CATCH_WEIGHT", "ARMIS_CATCH_WEIGHT", "OFF_LANDINGS_WEIGHT", "LANDINGS_MULT", "OFF_LANDINGS_VALUE")
	if(checkNms(object, rnms)==FALSE) stop("Check table \"Cl\" columns' size and names.")
	return(TRUE)
})

# CE
setGeneric("checkCEnms", function(object, ...){
	standardGeneric("checkCEnms")
	}
)

setMethod("checkCEnms", signature(object="data.frame"), function(object, ...){
	nms <- names(object)
	rnms <- c("RECORD_TYPE", "VESSEL_FLAG_COUNTRY", "YEAR", "QUARTER", "MONTH", "AREA", "RECTANGLE", "SUB_RECTANGLE", "FISHING_ACTIVITY_NAT", "FISHING_ACTIVITY_EU_L5", "FISHING_ACTIVITY_EU_L6", "HARBOUR", "NUMBER_OF_TRIPS", "NB_OF_SETS_HAULS", "FISHING_TIME", "KW_DAYS", "GT_DAYS", "DAYS_AT_SEA")
	if(checkNms(object, rnms)==FALSE) stop("Check table \"CE\" columns' size and names.")
	return(TRUE)
})

#====================================================================
# A set of methods to get and check PK in tables 
#====================================================================

# TR
setGeneric("checkTRpk", function(object, ...){
	standardGeneric("checkTRpk")
	}
)

setMethod("checkTRpk", signature(object="data.frame"), function(object, ...){
	identical(object[,1:6], unique(object[,1:6]))
})

# HH
setGeneric("checkHHpk", function(object, ...){
	standardGeneric("checkHHpk")
	}
)

setMethod("checkHHpk", signature(object="data.frame"), function(object, ...){
	identical(object[,1:7], unique(object[,1:7]))
})

# SL
setGeneric("checkSLpk", function(object, ...){
	standardGeneric("checkSLpk")
	}
)

setMethod("checkSLpk", signature(object="data.frame"), function(object, ...){
	identical(object[,1:14], unique(object[,1:14]))
})

# HL
setGeneric("checkHLpk", function(object, ...){
	standardGeneric("checkHLpk")
	}
)

setMethod("checkHLpk", signature(object="data.frame"), function(object, ...){
	identical(object[,1:15], unique(object[,1:15]))
})

# CA
setGeneric("checkCApk", function(object, ...){
	standardGeneric("checkCApk")
	}
)

setMethod("checkCApk", signature(object="data.frame"), function(object, ...){
	identical(object[,c(1:22)], unique(object[,c(1:22)]))
})

# CL
setGeneric("checkCLpk", function(object, ...){
	standardGeneric("checkCLpk")
	}
)

setMethod("checkCLpk", signature(object="data.frame"), function(object, ...){
	identical(object[,c(1:15)], unique(object[,c(1:15)]))
})

# CE
setGeneric("checkCEpk", function(object, ...){
	standardGeneric("checkCEpk")
	}
)

setMethod("checkCEpk", signature(object="data.frame"), function(object, ...){
	identical(object[,c(1:10)], unique(object[,c(1:10)]))
})

#====================================================================
# Methods to check data integrity 
#====================================================================

# data integrity 
setGeneric("checkDataIntegrity", function(target, current, ...){
	standardGeneric("checkDataIntegrity")
	}
)

setMethod("checkDataIntegrity", signature(target="data.frame", current="data.frame"), function(target, current, report=FALSE, ...){

	trg <- apply(target, 1, paste, collapse="")
	trg <- gsub("[[:space:]]","",trg)                                             #added MM 21/07/08
	current <- unique(current)
	if(sum(is.na(current))==ncol(current)) return(TRUE)
	crr <- apply(current, 1, paste, collapse="")
	crr <- gsub("[[:space:]]","",crr)                                             #added MM 21/07/08
	if(report==TRUE){
		current[!(crr %in% trg),]	
	} else {
		sum(crr %in% trg)==length(crr)
	}
})

#====================================================================
# check types of columns in tables 
#====================================================================

setGeneric("checkTys", function(object, tys, ...){
	standardGeneric("checkTys")
	}
)

setMethod("checkTys", signature("data.frame", "list"), function(object, tys, ...){
	n <- ncol(object)
	lst <- split(1:n, 1:n)
	lst <- lapply(lst, function(x) is(object[,x], tys[[x]]))
	identical(sum(unlist(lst)),n) 
})

