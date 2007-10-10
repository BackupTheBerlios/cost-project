#===================================
#
# EJ, 24/09/2007
# clData-class
#
#===================================

#====================================================================
# Class definition and validity check
#====================================================================

valclData <- function(object){

	cl <- object@cl

	# I will rely o the prototype to check col names and size. I'm not sure it's a good strategy !
	obj <- new("clData")
	cl0 <- obj@cl
	
	# check columns
	if(checkNms(cl, names(cl0))==FALSE) stop("Check slot candidate \"cl\" columns' size and names.")
	
	# check PK
	if(checkCLpk(cl)==FALSE) stop("Primary key not unique in slot candidate \"cl\".")

	# Everything is fine
	return(TRUE)
}

setClass("clData",
	representation(
		desc="character",
		cl="data.frame"
	),
	prototype(
		desc="my stock",
		cl=data.frame(
			landCtry=NA, # PK
			vslFlgCtry=NA, # PK
			year=NA, # PK
			quarter=NA, # PK 
			month=NA, # PK
			area=NA, # PK
			rectangle=NA, # PK 
			spp=NA, # PK 
			landCat=NA, # PK 
			commCatScl=NA, # PK
			commCat=NA, # PK
			foCatNat=NA, # PK
			foCatEu5=NA, # PK
			foCatEu6=NA, # PK
			unallocCatchWt=NA,
			misRepCatchWt=NA,
			landWt=NA,
			landMult=NA,
			landValue=NA)		
	),
	validity=valclData
)

#====================================================================
# Class constructor
#====================================================================
setGeneric("clData", function(cl, ...){
	standardGeneric("clData")
	}
)

setMethod("clData", signature("data.frame"), function(cl, desc="Unknown stock", ...){
	# create object and name columns properly 
	obj <- new("clData")
	names(cl) <- names(obj@cl)
	new("clData", cl=cl, desc=desc)
})

setMethod("clData", signature("missing"), function(desc="Unknown stock", ...){
	new("clData", desc=desc)
})

#====================================================================
# IO constructor
#====================================================================

setMethod("clData", signature("character"), function(cl, desc="Unknown stock", ...){

	# read CSV files
	# ToDo
	cl <- read.csv(cl)

	# check names are correct
	checkCLnms(cl)

	# remove record type 
	cl <- cl[,-1]

	# create object and name columns properly 
	obj <- new("clData")
	names(cl) <- names(obj@cl)
	new("clData", cl=cl, desc=desc)
})

#====================================================================
# Accessor functions
#====================================================================

setGeneric("cl", function(object, ...){
	standardGeneric("cl")
	}
)

setMethod("cl", signature("clData"), function(object, ...){
	object@cl
	}
)

setMethod("desc", signature("clData"), function(object, ...){
	object@desc
	}
)

#====================================================================
# 'Head' and 'Tail' functions  (setGeneric methods in 'csData-class.R')
#====================================================================


setMethod("head", signature("clData"), function(x, ...){
  object <- new("clData",desc=x@desc)
  object@cl <- head(x@cl)
  return(object)  
	}
)


setMethod("tail", signature("clData"), function(x, ...){
  object <- new("clData",desc=x@desc)
  object@cl <- tail(x@cl)
  return(object)  
	}
)

#====================================================================
# 'summary' function
#====================================================================

setMethod("summary", signature("clData"), function(object, ...){
  ll <- list()
  ll$desc <- object@desc
  ll$cl <- summary(object@cl)
  return(ll)  
	}
)

#====================================================================
# 'dim' function
#====================================================================

setMethod("dim", signature("clData"), function(x){
return(dim(x@cl))  
})

#====================================================================
# 'is.' function
#====================================================================

setGeneric("is.clData", function(object){
standardGeneric("is.clData")
})


setMethod("is.clData","ANY", function(object){
return(is(object)[1]=="clData")
})

