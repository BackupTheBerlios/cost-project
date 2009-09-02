#-------------------------------------------------------------------------------
# dbeOutputSim class and creators.
# Dorleta Garcia :: Azti-Tecnalia
# 29/04/2009 15:50:50
#-------------------------------------------------------------------------------


#library(COSTcore)
#
#====================================================================
# dbeObj : outcome object for all COSTdbe package methods ouputs
#====================================================================


setClass("OutputSim",
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
        desc        = "ObjectSim",
		species     = as.character(NA),
		catchCat    = as.character(NA),
		param       = as.character(NA),
		strataDesc  = strIni(),
		methodDesc  = as.character(NA),
		nSamp       = list(
                        len=data.frame(
                            sample = as.numeric(NA),
                            time=as.character(NA),
                            space=as.character(NA),
                            technical=as.character(NA),
                            value=as.numeric(NA)),
                        age=data.frame(
                            sample = as.numeric(NA),
                            time=as.character(NA),
                            space=as.character(NA),
                            value=as.numeric(NA))),
	                   nMeas       =list(      
                            len=data.frame(
                                sample = as.numeric(NA),
                                time=as.character(NA),
                                space=as.character(NA),
                                technical=as.character(NA),
                                value=as.numeric(NA)),
                            age=data.frame(
                                sample = as.numeric(NA),
                                time=as.character(NA),
                                space=as.character(NA),
                                value=as.numeric(NA))),
		              lenStruc=list(
                            estim=data.frame(    
                                sample = as.numeric(NA),  
                                time=as.character(NA),
                                space=as.character(NA),
                                technical=as.character(NA),
                                length=as.character(NA),
                                value=as.numeric(NA)),
                            rep=data.frame(
                                sample = as.numeric(NA),
                                time=as.character(NA),
                                space=as.character(NA),
                                technical=as.character(NA),
                                length=as.character(NA),
                                value=as.numeric(NA),
                                iter=as.numeric(NA))
                    ),
                    lenVar=data.frame(
                                sample = as.numeric(NA),      
                                time=as.character(NA),
                                space=as.character(NA),
                                technical=as.character(NA),
                                length=as.character(NA),
                                value=as.numeric(NA)),
                    lenNum=list(
                                
                                ci=data.frame(
                                sample = as.numeric(NA),
                                    time=as.character(NA),
                                    space=as.character(NA),
                                    technical=as.character(NA),
                                    length=as.character(NA),
                                    value=as.numeric(NA),
                                    inf=as.numeric(NA),
                                    sup=as.numeric(NA)),
                                cv=data.frame(
                                    sample = as.numeric(NA),
                                    time=as.character(NA),
                                    space=as.character(NA),
                                    technical=as.character(NA),
                                    length=as.character(NA),
                                    value=as.numeric(NA)),
                                DCRcvIndicator=data.frame(sample = as.numeric(NA),
                                                          value = as.numeric(NA))),    
	           ageStruc=list(
                        estim=data.frame( 
                            sample = as.numeric(NA),     
                            time=as.character(NA),
                            space=as.character(NA),
                            technical=as.character(NA),
                            age=as.character(NA),
                            value=as.numeric(NA)),
                        rep=data.frame(
                            sample = as.numeric(NA),
                            time=as.character(NA),
                            space=as.character(NA),
                            technical=as.character(NA),
                            age=as.character(NA),
                            value=as.numeric(NA),
                            iter=as.numeric(NA))),
		      ageVar=data.frame(
                        sample = as.numeric(NA),     
                        time=as.character(NA),
                        space=as.character(NA),
                        technical=as.character(NA),
                        age=as.character(NA),
                        value=as.numeric(NA)),
             ageNum=list(
                    ci=data.frame(
                        sample = as.numeric(NA),
                        time=as.character(NA),
                        space=as.character(NA),
                        technical=as.character(NA),
                        age=as.character(NA),
                        value=as.numeric(NA),
                        inf=as.numeric(NA),
                        sup=as.numeric(NA)),
                    cv=data.frame(
                        sample = as.numeric(NA),
                        time=as.character(NA),
                        space=as.character(NA),
                        technical=as.character(NA),
                        age=as.character(NA),
                        value=as.numeric(NA)),
                    DCRcvIndicator=data.frame(sample = as.numeric(NA),
                                                          value = as.numeric(NA))),  
            totalN=list(
                    estim=data.frame(
                        sample = as.numeric(NA),      
                        time=as.character(NA),
                        space=as.character(NA),
                        technical=as.character(NA),
                        value=as.numeric(NA)),
                    rep=data.frame(
                        sample = as.numeric(NA), 
                        time=as.character(NA),
                        space=as.character(NA),
                        technical=as.character(NA),
                        value=as.numeric(NA),
                        iter=as.numeric(NA))),
		  totalNvar=data.frame(
                    sample = as.numeric(NA),       
                    time=as.character(NA),
                    space=as.character(NA),
                    technical=as.character(NA),
                    value=as.numeric(NA)),
          totalNnum=list(
                    ci=data.frame(
                        sample = as.numeric(NA), 
                        time=as.character(NA),
                        space=as.character(NA),
                        technical=as.character(NA),
                        value=as.numeric(NA),
                        inf=as.numeric(NA),
                        sup=as.numeric(NA)),
                    cv=data.frame(
                        sample = as.numeric(NA), 
                        time=as.character(NA),
                        space=as.character(NA),
                        technical=as.character(NA),
                        value=as.numeric(NA)),
            DCRcvIndicator = data.frame(sample = as.numeric(NA),
                                                          value = as.numeric(NA))),    
		totalW=list(
            estim=data.frame(
                sample = as.numeric(NA),       
                time=as.character(NA),
                space=as.character(NA),
                technical=as.character(NA),
                value=as.numeric(NA)),
        rep=data.frame(
                sample = as.numeric(NA), 
                time=as.character(NA),
                space=as.character(NA),
                technical=as.character(NA),
                value=as.numeric(NA),
                iter=as.numeric(NA))),
        totalWvar=data.frame(
                sample = as.numeric(NA),       
                time=as.character(NA),
                space=as.character(NA),
                technical=as.character(NA),
                value=as.numeric(NA)),
        totalWnum=list(
            ci=data.frame(
                sample = as.numeric(NA), 
                time=as.character(NA),
                space=as.character(NA),
                technical=as.character(NA),
                value=as.numeric(NA),
                inf=as.numeric(NA),
                sup=as.numeric(NA)),
            cv=data.frame(
                sample = as.numeric(NA), 
                time=as.character(NA),
                space=as.character(NA),
                technical=as.character(NA),
                value=as.numeric(NA)),
        DCRcvIndicator=data.frame(sample = as.numeric(NA),
                                                          value = as.numeric(NA)))    
	)
)



#====================================================================
# 'dbeOutput' object constructor (initialization)
#====================================================================

ObjectSim <- function(desc, species, catchCat, param, strataDesc, methodDesc){

    if (missing(desc)) desc <- as.character(NA)
    if (missing(species)|all(is.na(species))) stop("Missing 'species' parameter!!")
    if (missing(catchCat)) catchCat <- "LAN"
    if ((length(catchCat)>1)|(all(is.na(catchCat)))) catchCat <- "LAN"
    if (missing(param)) param <- as.character(NA)
    if (missing(strataDesc)) strataDesc <- strIni()
    if (missing(methodDesc)) methodDesc <- as.character(NA)

new("OutputSim",desc=desc,species=species,catchCat=toupper(catchCat),param=param,strataDesc=strataDesc,methodDesc=methodDesc)
 	  }


#setGeneric("dbeOutputSim", function(obj, ...){
#	standardGeneric("dbeOutputSim")
#})
#setMethod("dbeOutputSim", signature(obj = 'missing'),   function(obj, desc, species, catchCat, param, strataDesc, methodDesc){ 
#    res <- dbeOutputSimt(desc=desc,species=species,catchCat=catchCat,param=param,strataDesc=strataDesc,methodDesc=methodDesc)
#return(res)})


