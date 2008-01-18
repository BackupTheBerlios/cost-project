#===================================
#
# EJ, 27/12/2007
# utils: assorted utility functions
#
#===================================

#====================================================================
# coerceDataFrameColumns
#====================================================================
setGeneric("coerceDataFrameColumns", function(object, refObject, ...){
	standardGeneric("coerceDataFrameColumns")
	}
)

setMethod("coerceDataFrameColumns", signature("data.frame", "data.frame"), function(object, refObject, ...){

	if(ncol(object)!=ncol(refObject)) stop("Both objects must have the same number of columns.\n")

	n <- ncol(object)
	for(i in 1:n){
		cls <- class(refObject[,i])
		if(cls=="factor"){
			v <- as(object[,i], "character")
			v <- as(object[,i], "factor")
		} else {
			v <- as(object[,i], cls)
		}
		object[,i] <- v
	}
	object	
})

