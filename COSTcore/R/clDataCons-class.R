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
			landCtry=as.factor(NA), # PK
			vslFlgCtry=as.factor(NA), # PK
			time=as.factor(NA), # PK
			space=as.factor(NA), # PK 
			technical=as.factor(NA), # PK
			spp=as.factor(NA), # PK 
			landCat=as.character(NA), # PK 
			commCatScl=as.character(NA), # PK
			commCat=as.character(NA), # PK
			unallocCatchWt=as.numeric(NA),
			misRepCatchWt=as.numeric(NA),
			landWt=as.numeric(NA),
			landMult=as.numeric(NA),
			landValue=as.numeric(NA))		
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

#====================================================================
# 'Head' and 'Tail' functions  (setGeneric methods in 'csData-class.R')
#====================================================================


setMethod("head", signature("clDataCons"), function(x, ...){
  object <- new("clDataCons",desc=x@desc)
  object@cl <- head(x@cl)
  return(object)  
	}
)


setMethod("tail", signature("clDataCons"), function(x, ...){
  object <- new("clDataCons",desc=x@desc)
  object@cl <- tail(x@cl)
  return(object)  
	}
)

#====================================================================
# 'summary' function
#====================================================================

setMethod("summary", signature("clDataCons"), function(object, ...){
  ll <- list()
  ll$desc <- object@desc
  ll$cl <- summary(object@cl)
  return(ll)  
	}
)

#====================================================================
# 'dim' function
#====================================================================

setMethod("dim", signature("clDataCons"), function(x){
	return(dim(x@cl))  
})

#====================================================================
# 'is.' function
#====================================================================

setGeneric("is.clDataCons", function(object){
	standardGeneric("is.clDataCons")
})


setMethod("is.clDataCons","ANY", function(object){
	return(is(object)[1]=="clDataCons")
})

#====================================================================
# rbind
#====================================================================

setMethod("rbind2", signature(x="clDataCons", y="clDataCons"), function(x,y){
	df0 <- rbind2(cl(x),cl(y))
	new("clDataCons", df0)
})

#====================================================================
# subset
#====================================================================

setMethod("subset", signature(x="clDataCons"), function(x,subset,...){
	e <- substitute(subset)
	df0 <- cl(x)	
	r <- eval(e, df0, parent.frame())
	new("clDataCons", df0[r,])
})

