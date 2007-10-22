#===================================
#
# EJ, 24/09/2007
# ceData-class
#
#===================================

#====================================================================
# Class definition and validity check
#====================================================================

valcecData <- function(object){

	ce <- object@ce

	# I will rely o the prototype to check col names and size. I'm not sure it's a good strategy !
	obj <- new("ceDataCons")
	ce0 <- obj@ce
	
	# check columns
	if(checkNms(ce, names(ce0))==FALSE) stop("Check slot candidate \"ce\" columns' size and names.")
	
	# check PK (ToDo)

	# Everything is fine
	return(TRUE)
}

setClass("ceDataCons",
	representation(
		desc="character",
		ce="data.frame"
	),
	prototype(
		desc="my stock",
		ce=data.frame(
			vslFlgCtry=NA, # PK
			time=NA, # PK
			space=NA, # PK 
			technical=NA, # PK
			trpNum=NA,
			foNum=NA,
			foDur=NA,
			effKwDays=NA,
			effGtDays=NA,
			daysAtSea=NA)		
	),
	validity=valcecData
)

#====================================================================
# Class constructor
#====================================================================
setGeneric("ceDataCons", function(object, ...){
	standardGeneric("ceDataCons")
	}
)

setMethod("ceDataCons", signature("ceDataVal"), function(object, ...){
	# create object and name columns properly 
	stop("Not implemented yet !")
})

setMethod("ceDataCons", signature("missing"), function(desc="Unknown stock", ...){
	new("ceDataCons", desc=desc)
})


