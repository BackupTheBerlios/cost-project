
#library(COSTcore)
#
#====================================================================
# dbeObj : outcome object for all COSTdbe package methods ouputs
#====================================================================


setClass("dbeOutputSim",
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
        desc        = "dbeObjectSim",
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

dbeObjectSim <- function(desc, species, catchCat, param, strataDesc, methodDesc, ...){

    if (missing(desc)) desc <- as.character(NA)
    if (missing(species)|all(is.na(species))) stop("Missing 'species' parameter!!")
    if (missing(catchCat)) catchCat <- "LAN"
    if ((length(catchCat)>1)|(all(is.na(catchCat)))) catchCat <- "LAN"
    if (missing(param)) param <- as.character(NA)
    if (missing(strataDesc)) strataDesc <- strIni()
    if (missing(methodDesc)) methodDesc <- as.character(NA)

new("dbeOutputSim",desc=desc,species=species,catchCat=toupper(catchCat),param=param,strataDesc=strataDesc,methodDesc=methodDesc)
 	  }
 	  
 	  
dbeSim2dbe <- function(dbeObjectSim, samples = unique(dbeObjectSim@nSamp$len$sample)){
      
      nn <- slotNames(dbeObject(species = '1'))
      bs  <- nn[which(!(sapply(nn, function(x) class(slot(dbeObject(species = '1'), x))) %in% c('list', 'data.frame')))]
      df  <- nn[which(sapply(nn, function(x) class(slot(dbeObject(species = '1'), x)))== 'data.frame')]
      lst <- nn[which(sapply(nn, function(x) class(slot(dbeObject(species = '1'), x)))== 'list')] 
      
      dbeObject <- dbeObject(species = dbeObjectSim@species)
       
      for(s in bs)  slot(dbeObject,s) <- slot(dbeObjectSim,s)
        
      for(s in df) slot(dbeObject,s) <- subset(slot(dbeObjectSim,s), sample %in% samples)[,-1]
      

      lst.df <- sapply(lst, function(x) lapply(slot(dbeObject(species = '1'),x), class))
      
      for(s in lst){
        for(sl in names(lst.df[[s]])){   
            slot(dbeObject,s)[[sl]] <- subset(slot(dbeObjectSim,s)[[sl]], sample %in% samples)[,-1]
        }
      }
      
      return(dbeObject)
}    

dbeSimNULL <- function(dbeObjectSim){
      
      df  <- names(which(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x)))== 'data.frame'))
      lst <- names(which(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x)))== 'list')) 
      lst.df <- sapply(lst, function(x) lapply(slot(dbeObject(species = '1'),x), class))

      res <- dbeObjectSim
              
      for(s in df) slot(res,s) <- slot(dbeObjectSim,s)[NULL,]
      
      for(s in lst){
        for(sl in names(lst.df[[s]])){   
            slot(res,s)[[sl]] <- slot(dbeObjectSim,s)[[sl]][NULL,]
        }
      }
      
      return(res)
}    
