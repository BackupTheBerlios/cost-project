#===================================
#
# EJ, 24/09/2007
# ceData-class
#
#===================================

#====================================================================
# Class definition and validity check
#====================================================================

valceData <- function(object){

	ce <- object@ce

	# I will rely o the prototype to check col names and size. I'm not sure it's a good strategy !
	obj <- new("ceData")
	ce0 <- obj@ce
	
	# check columns
	if(checkNms(ce, names(ce))==FALSE) stop("Check slot candidate \"ce\" columns' size and names.")
	
	# check PK
	if(checkCEpk(ce)==FALSE) stop("Primary key not unique in slot candidate \"ce\".")

	# Everything is fine
	return(TRUE)
}

setClass("ceData",
	representation(
		desc="character",
		ce="data.frame"
	),
	prototype(
		desc="my stock",
		ce=data.frame(
			vslFlgCtry=NA, # PK
			year=NA, # PK
			quarter=NA, # PK 
			month=NA, # PK
			area=NA, # PK
			rectangle=NA, # PK 
			foCatNat=NA, # PK
			foCatEu5=NA, # PK
			foCatEu6=NA, # PK
			trpNum=NA,
			foNum=NA,
			foDur=NA,
			effKwDays=NA,
			effGtDays=NA,
			daysAtSea=NA)		
	),
	validity=valceData
)

#====================================================================
# Class constructor
#====================================================================
setGeneric("ceData", function(ce, ...){
	standardGeneric("ceData")
	}
)

setMethod("ceData", signature("data.frame"), function(ce, desc="Unknown stock", ...){
	# create object and name columns properly 
	obj <- new("ceData")
	names(ce) <- names(obj@ce)
	new("ceData", ce=ce, desc=desc)
})

setMethod("ceData", signature("missing"), function(desc="Unknown stock", ...){
	new("ceData", desc=desc)
})

#====================================================================
# IO constructor
#====================================================================

setMethod("ceData", signature("character"), function(ce, desc="Unknown stock", ...){

	# read CSV files
	# ToDo
	ce <- read.csv(ce)

	# check names are correct
	checkCEnms(ce)

	# remove record type 
	ce <- ce[,-1]

	# create object and name columns properly 
	obj <- new("ceData")
	names(ce) <- names(obj@ce)
	new("ceData", ce=ce, desc=desc)
})

#====================================================================
# Accessor functions
#====================================================================

setGeneric("ce", function(object, ...){
	standardGeneric("ce")
	}
)

setMethod("ce", signature("ceData"), function(object, ...){
	object@ce
	}
)

setGeneric("desc", function(object, ...){
	standardGeneric("desc")
	}
)

setMethod("desc", signature("ceData"), function(object, ...){
	object@desc
	}
)
