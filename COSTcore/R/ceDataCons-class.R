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
			trpNum=as.numeric(NA),
			foNum=as.numeric(NA),
			foDur=as.numeric(NA),
			effKwDays=as.numeric(NA),
			effGtDays=as.numeric(NA),
			daysAtSea=as.numeric(NA)
		)),
	validity=valcecData
)




setGeneric("coerceCons", function(object, refObject, ...){
	standardGeneric("coerceCons")
	}
)

setMethod("coerceCons", signature("data.frame", "data.frame"), function(object, refObject, ...){

	if(ncol(object)!=ncol(refObject)) stop("Both objects must have the same number of columns.\n")

	n <- ncol(object)
	for(i in 1:n){
		cls <- class(refObject[,i])
		v <- as.character(object[,i])
		eval(parse('',text=paste("v <- as.",cls,"(v)",sep="")))
		object[,i] <- v
	}
	object	
})



#############
# "Stratif" #   =stratification definition for Cons objects creation
#====================================================================
# Class definition
#====================================================================

setClassUnion("NLchar",c("character","NULL"))
setClass("StratIni",representation(tempStrata="NLchar",spaceStrata="NLchar",techStrata="NLchar",sorting="NLchar"),
	                 prototype(tempStrata=NULL,spaceStrata=NULL,techStrata=NULL,sorting=NULL))		

#====================================================================
# Class constructor
#====================================================================

StratIni <- function(tempStrata=NULL,spaceStrata=NULL,techStrata=NULL,sorting=NULL) {                #sorting="catchCat" ou "commCat" ou "subSampcat"
new("StratIni",tempStrata=tempStrata,spaceStrata=spaceStrata,techStrata=techStrata,sorting=sorting)
}






#Aggregation tool
SpeedAgreg2 <- function(X,BY,FUN,...){
FactCar <- sapply(BY,as.character)
val <- apply(FactCar,1,function(x) paste(x,collapse="::"))
valAg <- aggregate(X,list(val=val),FUN,...)
tab <- as.data.frame(matrix(unlist(strsplit(as.character(valAg$val),"::")),ncol=length(BY),byrow=TRUE))
tab.ag <- data.frame(tab,valAg[,-1])
namBY <- names(BY) ; namX <- names(X)
if (is.null(namBY)) namBY <- rep("",length(BY)) ; if (is.null(namX)) namX <- rep("",length(X))
namBY[namBY==""] <- paste("c.",1:sum(namBY==""),sep="") ; namX[namX==""] <- paste("v.",1:sum(namX==""),sep="")
names(tab.ag) <- c(namBY,namX)
return(tab.ag)}









##====================================================================
## Class constructor
##====================================================================
#setGeneric("ceDataCons", function(object, ...){
#	standardGeneric("ceDataCons")
#	}
#)
#
#setMethod("ceDataCons", signature("ceDataVal"), function(object, ...){
#
#	ce <- ce(object)
#
#	#------------------------------------------------------------------------------
#	# time
#	#------------------------------------------------------------------------------
#	ce$time <- paste(ce$year, paste("Q", ce$quarter, sep=""), sep=".")
#
#	#------------------------------------------------------------------------------
#	# tech
#	#------------------------------------------------------------------------------
#	ce$technical <- apply(ce[,c("foCatNat","foCatEu5","foCatEu6")], 1,paste, collapse=".") 
#	
#	#------------------------------------------------------------------------------
#	# space
#	#------------------------------------------------------------------------------
#	ce$space <- apply(ce[,c("area","rect")], 1,paste, collapse=".") 
#	
#	#------------------------------------------------------------------------------
#	# create csDataCons
#	#------------------------------------------------------------------------------
#	cec <- ceDataCons()
#	ce <- ce[,match(names(ce(cec)),names(ce))]
#	new("ceDataCons", ce=ce)
#})
#
#setMethod("ceDataCons", signature("missing"), function(desc="Unknown stock", ...){
#	new("ceDataCons", desc=desc)
#})


setGeneric("ceDataCons", function(object,objStrat,...){
	standardGeneric("ceDataCons")
	}
)
		

setMethod("ceDataCons", signature("ceDataVal","StratIni"), function(object,objStrat,desc="Unknown stock",
                                                         TPrec=NULL,SPrec=NULL,TCrec=NULL,...){  #ex: TPrec=list(from=c("1","2","3","4"),to=c("5","5","6","6"))

tempStrata <- objStrat@tempStrata ; spaceStrata <- objStrat@spaceStrata ; techStrata <- objStrat@techStrata
if (techStrata=="commCat") stop("effort object do not match with market category sampling strategy")
CE <- object@ce 
CE$semester <- ceiling(CE$quarter/2)      
if (is.null(tempStrata)) {CE$time <- NA ; TPrec <- NULL} else CE$time <- CE[,tempStrata]     
if (is.null(spaceStrata)) {CE$space <- NA ; SPrec <- NULL} else CE$space <- CE[,spaceStrata]
if (is.null(techStrata)) {CE$technical <- NA ; TCrec <- NULL} else CE$technical <- CE[,techStrata]

#on recode si besoin est
if (!is.null(TPrec)) {Typ <- class(CE$time) ; CE$time <- factor(CE$time) ; Lev <- levels(CE$time)[!levels(CE$time)%in%TPrec$from]
                      CE$time <- factor(CE$time,levels=c(Lev,TPrec$from),labels=c(Lev,TPrec$to)) ; eval(parse('',text=paste("CE$time <- as.",Typ,"(as.character(CE$time))",sep="")))}
if (!is.null(SPrec)) {Typ <- class(CE$space) ; CE$space <- factor(CE$space) ; Lev <- levels(CE$space)[!levels(CE$space)%in%SPrec$from]
                      CE$space <- factor(CE$space,levels=c(Lev,SPrec$from),labels=c(Lev,SPrec$to)) ; eval(parse('',text=paste("CE$space <- as.",Typ,"(as.character(CE$space))",sep="")))}
if (!is.null(TCrec)) {Typ <- class(CE$technical) ; CE$technical <- factor(CE$technical) ; Lev <- levels(CE$technical)[!levels(CE$technical)%in%TCrec$from]
                      CE$technical <- factor(CE$technical,levels=c(Lev,TCrec$from),labels=c(Lev,TCrec$to)) ; eval(parse('',text=paste("CE$technical <- as.",Typ,"(as.character(CE$technical))",sep="")))}
                        

csc <- new("ceDataCons")
	ce <- CE[,match(names(csc@ce),names(CE))] ; rownames(ce) <- 1:nrow(ce)  
new("ceDataCons",desc=desc,ce=coerceCons(ce,csc@ce))
	
})



setMethod("ceDataCons", signature("ceDataVal","missing"), function(object,desc="Unknown stock", ...){

	ceDataCons(object,StratIni(),desc=desc,...)
})

setMethod("ceDataCons", signature("missing","missing"), function(desc="Unknown stock", ...){

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

