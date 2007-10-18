#===================================
#
# EJ, 24/09/2007
# csData-class
#
#===================================

#====================================================================
# Class definition and validity check
#====================================================================

setClass("csDataVal", contains="csData")

#====================================================================
# Class constructor
#====================================================================
setGeneric("csDataVal", function(object, ...){
	standardGeneric("csDataVal")
	}
)

setMethod("csDataVal", signature("csData"), function(object, ...){
	new("csDataVal", tr=tr(object), hh=hh(object), sl=sl(object), hl=hl(object), ca=ca(object), desc=desc(object))
})

setMethod("csDataVal", signature("missing"), function(desc="Unknown stock", ...){
	new("csDataVal", desc=desc)
})


