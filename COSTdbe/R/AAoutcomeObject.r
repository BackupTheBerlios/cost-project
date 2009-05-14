
#library(COSTcore)
#
#====================================================================
# dbeObj : outcome object for all COSTdbe package methods ouputs
#====================================================================


setClass("dbeOutput",
	representation(
    desc="character",                        #descriptor
    species="character",                     #recall of SL$spp (+ SL$sex)
		catchCat="character",                    #recall of the catch category (discards/landings)
		param="character",                       #recall of the parameter estimated (N, W, maturity, sex-ratio,...)
		strataDesc="strIni",                     #time, space and technical stratification considered
		methodDesc="character",                  #recall of the method (analytical, bootstrap, bayesian)
		nSamp="list",                            #number of samples
		nMeas="list",                            #number of individual measures
		lenStruc="list",                         #estimates of the length structure (param-at-length)
		lenVar="data.frame",                     #estimates of the variance of '$lenStruc'
		lenNum="list",                           #further numerical data about length structure (ex: ci, cv) 
		ageStruc="list",                         #estimates of the age structure (param-at-age)
		ageVar="data.frame",                     #estimates of the variance of '$ageStruc'
		ageNum="list",                           #further numerical data about age structure (ex: ci, cv) 
		totalN="list",                           #estimates of the total number of the parameters
		totalNvar="data.frame",                  #estimates of the variance of '$totalN'
		totalNnum="list",                        #further numerical data about total numbers (ex: ci, cv) 
		totalW="list",                           #estimates of the total weight of the parameters
		totalWvar="data.frame",                  #estimates of the variance of '$totalW'
		totalWnum="list"                         #further numerical data about total weights (ex: ci, cv) 
	),
	prototype(
    desc="dbeObject",
		species=as.character(NA),
		catchCat=as.character(NA),
		param=as.character(NA),
		strataDesc=strIni(),
		methodDesc=as.character(NA),
		nSamp=list(
      len=data.frame(
        time=as.character(NA),
        space=as.character(NA),
        technical=as.character(NA),
        value=as.numeric(NA)),
      age=data.frame(
        time=as.character(NA),
        space=as.character(NA),
        value=as.numeric(NA))
      ),
		nMeas=list(      
      len=data.frame(
        time=as.character(NA),
        space=as.character(NA),
        technical=as.character(NA),
        value=as.numeric(NA)),
      age=data.frame(
        time=as.character(NA),
        space=as.character(NA),
        value=as.numeric(NA))
      ),
		lenStruc=list(
        estim=data.frame(      
          time=as.character(NA),
          space=as.character(NA),
          technical=as.character(NA),
          length=as.character(NA),
          value=as.numeric(NA)),
        rep=data.frame(
          time=as.character(NA),
          space=as.character(NA),
          technical=as.character(NA),
          length=as.character(NA),
          value=as.numeric(NA),
          iter=as.numeric(NA))
          ),
		lenVar=data.frame(      
      time=as.character(NA),
      space=as.character(NA),
      technical=as.character(NA),
      length=as.character(NA),
      value=as.numeric(NA)),
    lenNum=list(
        ci=data.frame(
          time=as.character(NA),
          space=as.character(NA),
          technical=as.character(NA),
          length=as.character(NA),
          value=as.numeric(NA),
          inf=as.numeric(NA),
          sup=as.numeric(NA)),
        cv=data.frame(
          time=as.character(NA),
          space=as.character(NA),
          technical=as.character(NA),
          length=as.character(NA),
          value=as.numeric(NA)),
        DCRcvIndicator=as.numeric(NA)
        ),    
		ageStruc=list(
        estim=data.frame(      
          time=as.character(NA),
          space=as.character(NA),
          technical=as.character(NA),
          age=as.character(NA),
          value=as.numeric(NA)),
        rep=data.frame(
          time=as.character(NA),
          space=as.character(NA),
          technical=as.character(NA),
          age=as.character(NA),
          value=as.numeric(NA),
          iter=as.numeric(NA))
          ),
		ageVar=data.frame(      
      time=as.character(NA),
      space=as.character(NA),
      technical=as.character(NA),
      age=as.character(NA),
      value=as.numeric(NA)),
    ageNum=list(
        ci=data.frame(
          time=as.character(NA),
          space=as.character(NA),
          technical=as.character(NA),
          age=as.character(NA),
          value=as.numeric(NA),
          inf=as.numeric(NA),
          sup=as.numeric(NA)),
        cv=data.frame(
          time=as.character(NA),
          space=as.character(NA),
          technical=as.character(NA),
          age=as.character(NA),
          value=as.numeric(NA)),
        DCRcvIndicator=as.numeric(NA)
        ),    
		totalN=list(
        estim=data.frame(      
          time=as.character(NA),
          space=as.character(NA),
          technical=as.character(NA),
          value=as.numeric(NA)),
        rep=data.frame(
          time=as.character(NA),
          space=as.character(NA),
          technical=as.character(NA),
          value=as.numeric(NA),
          iter=as.numeric(NA))
          ),
		totalNvar=data.frame(      
      time=as.character(NA),
      space=as.character(NA),
      technical=as.character(NA),
      value=as.numeric(NA)),
    totalNnum=list(
        ci=data.frame(
          time=as.character(NA),
          space=as.character(NA),
          technical=as.character(NA),
          value=as.numeric(NA),
          inf=as.numeric(NA),
          sup=as.numeric(NA)),
        cv=data.frame(
          time=as.character(NA),
          space=as.character(NA),
          technical=as.character(NA),
          value=as.numeric(NA)),
        DCRcvIndicator=as.numeric(NA)
        ),    
		totalW=list(
        estim=data.frame(      
          time=as.character(NA),
          space=as.character(NA),
          technical=as.character(NA),
          value=as.numeric(NA)),
        rep=data.frame(
          time=as.character(NA),
          space=as.character(NA),
          technical=as.character(NA),
          value=as.numeric(NA),
          iter=as.numeric(NA))
          ),
		totalWvar=data.frame(      
      time=as.character(NA),
      space=as.character(NA),
      technical=as.character(NA),
      value=as.numeric(NA)),
    totalWnum=list(
        ci=data.frame(
          time=as.character(NA),
          space=as.character(NA),
          technical=as.character(NA),
          value=as.numeric(NA),
          inf=as.numeric(NA),
          sup=as.numeric(NA)),
        cv=data.frame(
          time=as.character(NA),
          space=as.character(NA),
          technical=as.character(NA),
          value=as.numeric(NA)),
        DCRcvIndicator=as.numeric(NA)
        )    
	)
)



#====================================================================
# 'dbeOutput' object constructor (initialization)
#====================================================================

dbeObject <- function(desc, species, catchCat, param, strataDesc, methodDesc, ...){

if (missing(desc)) desc <- as.character(NA)
if (missing(species)|all(is.na(species))) stop("Missing 'species' parameter!!")
if (missing(catchCat)) catchCat <- "LAN"
if ((length(catchCat)>1)|(all(is.na(catchCat)))) catchCat <- "LAN"
if (missing(param)) param <- as.character(NA)
if (missing(strataDesc)) strataDesc <- strIni()
if (missing(methodDesc)) methodDesc <- as.character(NA)

new("dbeOutput",desc=desc,species=species,catchCat=toupper(catchCat),param=param,strataDesc=strataDesc,methodDesc=methodDesc)
 	  }
 	  

#====================================================================
# 'rbind2' method for 'dbeOutput' object 
#====================================================================

 
setMethod("rbind2", signature(x="dbeOutput", y="dbeOutput"), function(x,y){

#subfunction
rBind2 <- function(tab1,tab2){
if (all(is.na(tab1))) {
  res <- tab2
} else {
  if (all(is.na(tab2))) {
    res <- tab1
  } else {
    res <- rbind2(tab1,tab2)
    rownames(res) <- 1:nrow(res)
  }
}
return(res)
}

elt <- c("nSamp$len","nSamp$age","nMeas$len","nMeas$age",
         "lenStruc$estim","lenStruc$rep","lenVar","lenNum$ci","lenNum$cv","lenNum$DCRcvIndicator",
         "ageStruc$estim","ageStruc$rep","ageVar","ageNum$ci","ageNum$cv","ageNum$DCRcvIndicator",
         "totalN$estim","totalN$rep","totalNvar","totalNnum$ci","totalNnum$cv","totalNnum$DCRcvIndicator",
         "totalW$estim","totalW$rep","totalWvar","totalWnum$ci","totalWnum$cv","totalWnum$DCRcvIndicator")

invisible(sapply(elt,function(z) eval(parse('',text=paste("x@",z," <<- rBind2(x@",z,",y@",z,")",sep=""))))) 
#'DCRcvIndicator' elements are set to NA (could be recalculated from x & y : IND = [ INDx * sum(ESTx) + INDy * sum(ESTy) ] / sum (ESTx + ESTy) <<- left to be implemented)
x@lenNum$DCRcvIndicator <- x@ageNum$DCRcvIndicator <- x@totalNnum$DCRcvIndicator <- x@totalWnum$DCRcvIndicator <- NA

return(x)
})


                   
#====================================================================
# '+' method for 'dbeOutput' object 
#====================================================================

	
   
     
setMethod("+", signature(e1="dbeOutput", e2="dbeOutput"), function(e1,e2){

#subfunction
addDBE <- function(tab1,tab2) {                
if (all(is.na(tab1))) {
  res <- tab2
} else {
if (all(is.na(tab2))) {
  res <- tab1
} else {
  TAB <- rbind(tab1,tab2)
  res <- aggregate(TAB$value,as.list(TAB[,rev(names(TAB)[-match("value",names(TAB))])]),sum,na.rm=TRUE)
  names(res)[ncol(res)] <- "value"
}}
return(res[,names(tab1)])
} 


elt <- c("nSamp$len","nSamp$age","nMeas$len","nMeas$age",
         "lenStruc$estim","lenStruc$rep","lenVar",
         "ageStruc$estim","ageStruc$rep","ageVar",
         "totalN$estim","totalN$rep","totalNvar",
         "totalW$estim","totalW$rep","totalWvar")

invisible(sapply(elt,function(z) eval(parse('',text=paste("e1@",z," <<- addDBE(e1@",z,",e2@",z,")",sep=""))))) 

eltNA <- c("lenNum$ci","lenNum$cv","lenNum$DCRcvIndicator",
           "ageNum$ci","ageNum$cv","ageNum$DCRcvIndicator",
           "totalNnum$ci","totalNnum$cv","totalNnum$DCRcvIndicator",
           "totalWnum$ci","totalWnum$cv","totalWnum$DCRcvIndicator")

invisible(sapply(eltNA,function(z) eval(parse('',text=paste("e1@",z," <<- new(\"dbeOutput\")@",z,sep=""))))) 

return(e1)

})



