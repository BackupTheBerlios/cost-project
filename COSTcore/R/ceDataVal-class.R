#===================================
#
# EJ, 24/09/2007
# ceData-class
#
#===================================

#====================================================================
# Class definition and validity check
#====================================================================

valcevData <- function(object){

	ce <- object@ce

	# I will rely o the prototype to check col names and size. I'm not sure it's a good strategy !
	obj <- new("ceData")
	ce0 <- obj@ce
	
	# check columns
	if(checkNms(ce, names(ce0))==FALSE) stop("Check slot candidate \"ce\" columns' size and names.")
	
	# check PK
	if(checkCEpk(ce)==FALSE) stop("Primary key not unique in slot candidate \"ce\".")

	# Everything is fine
	return(TRUE)
}

setClass("ceDataVal", contains="ceData", validity=valcevData)

#====================================================================
# Class constructor
#====================================================================
setGeneric("ceDataVal", function(object, ...){
	standardGeneric("ceDataVal")
	}
)

setMethod("ceDataVal", signature("ceData"), function(object, ...){
	new("ceDataVal", ce=ce(object), desc=desc(object))
})

setMethod("ceDataVal", signature("missing"), function(desc="Unknown stock", ...){
	new("ceDataVal", desc=desc)
})


