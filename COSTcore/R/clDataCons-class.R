#===================================
#
# EJ, 24/09/2007
# clData-class
#
#===================================

#====================================================================
# Class definition and validity check
#====================================================================

valclcData <- function(object){

	cl <- object@cl

	# I will rely o the prototype to check col names and size. I'm not sure it's a good strategy !
	obj <- new("clDataCons")
	cl0 <- obj@cl
	
	# check columns
	if(checkNms(cl, names(cl0))==FALSE) stop("Check slot candidate \"ce\" columns' size and names.")
	
	# check PK (ToDo)

	# Everything is fine
	return(TRUE)
}

setClass("clDataCons", 
	representation(
		desc="character",
		cl="data.frame"
	),
	prototype(
		desc="my stock",
		cl=data.frame(
			landCtry=NA, # PK
			vslFlgCtry=NA, # PK
			time=NA, # PK
			space=NA, # PK 
			technical=NA, # PK
			spp=NA, # PK 
			landCat=NA, # PK 
			commCatScl=NA, # PK
			commCat=NA, # PK
			unallocCatchWt=NA,
			misRepCatchWt=NA,
			landWt=NA,
			landMult=NA,
			landValue=NA)		
	),
	validity=valclcData
)

#====================================================================
# Class constructor
#====================================================================
setGeneric("clDataCons", function(object, ...){
	standardGeneric("clDataCons")
	}
)

setMethod("clDataCons", signature("clDataVal"), function(object, ...){
	# create object and name columns properly 
	stop("Not implemented yet !")
})

setMethod("clDataCons", signature("missing"), function(desc="Unknown stock", ...){
	new("clDataCons", desc=desc)
})


