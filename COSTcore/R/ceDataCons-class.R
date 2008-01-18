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
			vslFlgCtry=as.factor(NA), # PK
			time=as.factor(NA), # PK
			space=as.factor(NA), # PK 
			technical=as.factor(NA), # PK
			trpNum=as.numeric(NA),
			foNum=as.numeric(NA),
			foDur=as.numeric(NA),
			effKwDays=as.numeric(NA),
			effGtDays=as.numeric(NA),
			daysAtSea=as.numeric(NA)
		)),
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

#====================================================================
# Accessor functions
#====================================================================

setMethod("ce", signature("ceDataCons"), function(object, ...){
	object@ce
	}
)

setMethod("desc", signature("ceDataCons"), function(object, ...){
	object@desc
	}
)

#====================================================================
# 'Head' and 'Tail' functions
#====================================================================


setMethod("head", signature("ceDataCons"), function(x, ...){
  object <- new("ceDataCons",desc=x@desc)
  object@ce <- head(x@ce)
  return(object)  
	}
)


setMethod("tail", signature("ceDataCons"), function(x, ...){
  object <- new("ceDataCons",desc=x@desc)
  object@ce <- tail(x@ce)
  return(object)  
	}
)

#====================================================================
# 'summary' function
#====================================================================

setMethod("summary", signature("ceDataCons"), function(object, ...){
  ll <- list()
  ll$desc <- object@desc
  ll$ce <- summary(object@ce)
  return(ll)  
	}
)


#====================================================================
# 'dim' function
#====================================================================

setMethod("dim", signature("ceDataCons"), function(x){
	return(dim(x@ce))  
})

#====================================================================
# 'is.' function
#====================================================================

setGeneric("is.ceDataCons", function(object){
	standardGeneric("is.ceDataCons")
})


setMethod("is.ceDataCons","ANY", function(object){
	return(is(object)[1]=="ceDataCons")
})

#====================================================================
# rbind
#====================================================================

setMethod("rbind2", signature(x="ceDataCons", y="ceDataCons"), function(x,y){
	df0 <- rbind2(ce(x),ce(y))
	new("ceDataCons", ce=df0)
})

#====================================================================
# subset
#====================================================================

setMethod("subset", signature(x="ceDataCons"), function(x,subset,...){
	e <- substitute(subset)
	df0 <- ce(x)	
	r <- eval(e, df0, parent.frame())
	new("ceDataCons", df0[r,])
})

