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
# 			year=as.numeric(NA), # PK	=> time
# 			quarter=as.numeric(NA), # PK	=> time 
# 			month=as.numeric(NA), # PK	=> time
# 			area=as.character(NA), # PK	=> space
# 			rect=as.character(NA), # PK 	=> space
# 			foCatNat=as.character(NA), # PK	=> tech
# 			foCatEu5=as.character(NA), # PK	=> tech
# 			foCatEu6=as.character(NA), # PK	=> tech
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
setGeneric("clDataCons", function(object,objStrat,...){
	standardGeneric("clDataCons")
	}
)

#setMethod("clDataCons", signature("clDataVal"), function(object, ...){
#
#	cl <- cl(object)
#
#	#------------------------------------------------------------------------------
#	# time
#	#------------------------------------------------------------------------------
#	cl$time <- paste(cl$year, paste("Q", cl$quarter, sep=""), sep=".")
#
#	#------------------------------------------------------------------------------
#	# tech
#	#------------------------------------------------------------------------------
#	cl$technical <- apply(cl[,c("foCatNat","foCatEu5","foCatEu6")], 1,paste, collapse=".") 
#	
#	#------------------------------------------------------------------------------
#	# space
#	#------------------------------------------------------------------------------
#	cl$space <- apply(cl[,c("area","rect")], 1,paste, collapse=".") 
#	
#	#------------------------------------------------------------------------------
#	# create csDataCons
#	#------------------------------------------------------------------------------
#	clc <- clDataCons()
#	cl <- cl[,match(names(cl(clc)),names(cl))]
#	new("clDataCons", cl=cl)
#})
#
#setMethod("clDataCons", signature("missing"), function(desc="Unknown stock", ...){
#	new("clDataCons", desc=desc)
#})


setMethod("clDataCons", signature("clDataVal","strIni"), function(object,
                                                                  objStrat,
                                                                  desc="Unknown stock",
                                                                  ...){  

timeStrata <- objStrat@timeStrata
spaceStrata <- objStrat@spaceStrata 
techStrata <- objStrat@techStrata
tpRec <- objStrat@tpRec
spRec <- objStrat@spRec
tcRec <- objStrat@tcRec

CL <- object@cl 
CL$semester <- ceiling(CL$quarter/2)      

#-------------------------------------------------------------------------------
# Creation of the 3 stratification fields in cl
#-------------------------------------------------------------------------------

if (is.na(timeStrata)) {
  CL$time <- NA
  tpRec <- as.list(NA)
} else {
  CL$time <- CL[,timeStrata]}    
   
if (is.na(spaceStrata)) {
  CL$space <- NA 
  spRec <- as.list(NA)
} else {
  CL$space <- CL[,spaceStrata]}

if (is.na(techStrata)) {
  CL$technical <- NA 
  tcRec <- as.list(NA)
} else {
  CL$technical <- CL[,techStrata]}


#-------------------------------------------------------------------------------
# Recoding the 3 new stratification fields following user post-stratification          <<<- there's surely a more simple way to do this
#-------------------------------------------------------------------------------

if (!is.na(tpRec[1])) {
  Typ <- class(CL$time)
  CL$time <- factor(CL$time)
  Lev <- levels(CL$time)[!levels(CL$time)%in%tpRec$from]
  CL$time <- factor(CL$time,levels=c(Lev,tpRec$from),labels=c(Lev,tpRec$to))
  eval(parse('',text=paste("CL$time <- as.",Typ,"(as.character(CL$time))",sep="")))}
  
if (!is.na(spRec[1])) {
  Typ <- class(CL$space)
  CL$space <- factor(CL$space)
  Lev <- levels(CL$space)[!levels(CL$space)%in%spRec$from]
  CL$space <- factor(CL$space,levels=c(Lev,spRec$from),labels=c(Lev,spRec$to))
  eval(parse('',text=paste("CL$space <- as.",Typ,"(as.character(CL$space))",sep="")))}
  
if (!is.na(tcRec[1])) {
  Typ <- class(CL$technical)
  CL$technical <- factor(CL$technical)
  Lev <- levels(CL$technical)[!levels(CL$technical)%in%tcRec$from]
  CL$technical <- factor(CL$technical,levels=c(Lev,tcRec$from),labels=c(Lev,tcRec$to))
  eval(parse('',text=paste("CL$technical <- as.",Typ,"(as.character(CL$technical))",sep="")))}
                        
#-------------------------------------------------------------------------------
# Finally, creation of 'clDataCons' object and use of 'coerceCons' function to convert columns
#-------------------------------------------------------------------------------

csc <- new("clDataCons")
cl <- CL[,match(names(csc@cl),names(CL))]
rownames(cl) <- 1:nrow(cl)  
new("clDataCons", desc=desc,cl=coerceCons(cl,csc@cl))
})
	





setMethod("clDataCons", signature("clDataVal","missing"), function(object,desc="Unknown stock", ...){

	clDataCons(object,strIni(),desc=desc,...)
})





setMethod("clDataCons", signature("missing","missing"), function(desc="Unknown stock", ...){

	new("clDataCons", desc=desc)
})	







#====================================================================
# Accessor functions
#====================================================================

setMethod("cl", signature("clDataCons"), function(object, ...){
	object@cl
	}
)

setMethod("desc", signature("clDataCons"), function(object, ...){
	object@desc
	}
)

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

