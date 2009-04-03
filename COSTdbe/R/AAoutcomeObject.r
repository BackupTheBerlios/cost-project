
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
 	  
 	  
