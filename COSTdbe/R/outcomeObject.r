library(COSTcore)

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
		nSamp="data.frame",                      #number of samples
		nMes="data.frame",                       #number of individual measures
		lenStruc="list",                         #estimates of the length structure (param-at-length)
		lenVar="data.frame",                     #estimates of the variance of '$lenStruc'
		ageStruc="list",                         #estimates of the age structure (param-at-age)
		ageVar="data.frame",                     #estimates of the variance of '$ageStruc'
		totalN="list",                           #estimates of the total number of the parameters
		totalNvar="data.frame",                  #estimates of the variance of '$totalN'
		totalW="list",                           #estimates of the total weight of the parameters
		totalWvar="data.frame"                   #estimates of the variance of '$totalW'
	),
	prototype(
    desc="dbeObject",
		species=as.character(NA),
		catchCat=as.character(NA),
		param=as.character(NA),
		strataDesc=strIni(),
		methodDesc=as.character(NA),
		nSamp=data.frame(
      time=as.character(NA),
      space=as.character(NA),
      technical=as.character(NA),
      value=as.numeric(NA)),
		nMes=data.frame(      
      time=as.character(NA),
      space=as.character(NA),
      technical=as.character(NA),
      value=as.numeric(NA)),
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
      value=as.numeric(NA))
	)
)



#====================================================================
# 'dbeOutput' object constructor (initialization)
#====================================================================

dbeObject <- function(desc, species, catchCat, param, strataDesc, methodDesc, ...){

if (missing(desc)) desc <- as.character(NA)
if (missing(species)|all(is.na(species))) stop("Missing 'species' parameter!!")
if (all(is.na(catchCat))) stop("Missing 'catchCat' parameter!!")
if (missing(catchCat)|(length(catchCat)>1)|(!all(tolower(catchCat)%in%c("dis","lan")))) stop("Wrong 'catchCat' parameter!!")
if (missing(param)) param <- as.character(NA)
if (missing(strataDesc)) strataDesc <- strIni()
if (missing(methodDesc)) methodDesc <- as.character(NA)

new("dbeOutput",desc=desc,species=species,catchCat=toupper(catchCat),param=param,strataDesc=strataDesc,methodDesc=methodDesc)
 	  }
 	  
 	  
