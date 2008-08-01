library(COSTcore)

#====================================================================
# dbeObj : outcome object for all COSTdbe package methods ouputs
#====================================================================

setClassUnion("numArray", c("numeric", "array"))                 #maybe not a good idea!!

setClass("dbeOutput",
	representation(
    desc="character",                      #descriptor
    species="character",                   #recall of SL$spp (+ SL$sex)
		catchCat="character",                  #recall of the catch category (discards/landings)
		param="character",                     #recall of the parameter estimated (N, W, maturity, sex-ratio,...)
		strataDesc="strIni",                   #time, space and technical stratification considered
		methodDesc="character",                #recall of the method (analytical, bootstrap, bayesian)
		nSamp="numArray",                      #number of samples
		nMes="numArray",                       #number of individual measures
		lenStruc="numArray",                   #estimates of the length structure (param-at-length)
		lenVar="numArray",                     #estimates of the variance of '$lenStruc'
		ageStruc="numArray",                   #estimates of the age structure (param-at-age)
		ageVar="numArray",                     #estimates of the variance of '$ageStruc'
		totalN="numArray",                     #estimates of the total number of the parameters
		totalNvar="numArray",                  #estimates of the variance of '$totalN'
		totalW="numArray",                     #estimates of the total weight of the parameters
		totalWvar="numArray"                   #estimates of the variance of '$totalW'
	),
	prototype(
    desc="dbeObject",
		species=as.character(NA),
		catchCat=as.character(NA),
		param=as.character(NA),
		strataDesc=strIni(),
		methodDesc=as.character(NA),
		nSamp=as.numeric(NA),
		nMes=as.numeric(NA),
		lenStruc=as.numeric(NA),
		lenVar=as.numeric(NA),
		ageStruc=as.numeric(NA),
		ageVar=as.numeric(NA),
		totalN=as.numeric(NA),
		totalNvar=as.numeric(NA),
		totalW=as.numeric(NA),
		totalWvar=as.numeric(NA)
	)
)






