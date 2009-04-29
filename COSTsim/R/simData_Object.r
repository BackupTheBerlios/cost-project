
#library(COSTcore)
#
#====================================================================
# dbeObj : outcome object for all COSTdbe package methods ouputs
#====================================================================


setClass("simData",
	representation(
        desc        = "character",                     # descriptor
        species     = "character",                     # recall of SL$spp (+ SL$sex)
        samples     = "list",                          # list with simulates 'costData' objects.
        initial.fit = "list",
		setup.args  = "list",                          # list with setup parameters
		burnin      = "numeric",                     #recall of the parameter estimated (N, W, maturity, sex-ratio,...)
		nmcmc       = "numeric",                        #time, space and technical stratification considered
		l.int       = "numeric",                     #recall of the method (analytical, bootstrap, bayesian)
		Int         = "list",                          #number of samples
		Slp         = "list",                          #number of individual measures
		landings    = "numeric",                          #estimates of the length structure (param-at-length)
		nHaul       = "integer",                    #estimates of the variance of '$lenStruc'
		nseas       = "numeric"                          #further numerical data about length structure (ex: ci, cv) )
		))
#	prototype(
#         desc        = "simDataObject",
#         samples     = list(),                        
#         initial.fit = NA,
#		 setup.data  = list(),                        
#		 burnin      = NA,                     
#		 nmcmc       = NA,                     
#		 l.int       = list(),               
#		 Int         = list(),                       
#		 Slp         = list(),                      
#		 landings    = NA,                   
#		 nHaul       = NA,                  
#		 nseas       = NA                    
#   ))


setClass("simDataVal",
	representation(
        desc        = "character",                     # descriptor
        species     = "character",                     # recall of SL$spp (+ SL$sex)
        samples     = "list",                          # list with simulates 'costData' objects.
        initial.fit = "list",
		setup.args  = "list",                          # list with setup parameters
		burnin      = "numeric",                     #recall of the parameter estimated (N, W, maturity, sex-ratio,...)
		nmcmc       = "numeric",                        #time, space and technical stratification considered
		l.int       = "numeric",                     #recall of the method (analytical, bootstrap, bayesian)
		Int         = "list",                          #number of samples
		Slp         = "list",                          #number of individual measures
		landings    = "numeric",                          #estimates of the length structure (param-at-length)
		nHaul       = "integer",                    #estimates of the variance of '$lenStruc'
		nseas       = "numeric"                          #further numerical data about length structure (ex: ci, cv) )
		))
		
setClass("simDataCons",
	representation(
        desc        = "character",                     # descriptor
        species     = "character",                     # recall of SL$spp (+ SL$sex)
        samples     = "list",                          # list with simulates 'costData' objects.
        initial.fit = "list",
		setup.args  = "list",                          # list with setup parameters
		burnin      = "numeric",                     #recall of the parameter estimated (N, W, maturity, sex-ratio,...)
		nmcmc       = "numeric",                        #time, space and technical stratification considered
		l.int       = "numeric",                     #recall of the method (analytical, bootstrap, bayesian)
		Int         = "list",                          #number of samples
		Slp         = "list",                          #number of individual measures
		landings    = "numeric",                          #estimates of the length structure (param-at-length)
		nHaul       = "integer",                    #estimates of the variance of '$lenStruc'
		nseas       = "numeric"                          #further numerical data about length structure (ex: ci, cv) )
		))
        		

#	prototype(
#         desc        = "simDataObject",
#         samples     = list(),                        
#         initial.fit = NA,
#		 setup.data  = list(),                        
#		 burnin      = NA,                     
#		 nmcmc       = NA,                     
#		 l.int       = list(),               
#		 Int         = list(),                       
#		 Slp         = list(),                      
#		 landings    = NA,                   
#		 nHaul       = NA,                  
#		 nseas       = NA                    
#   ))


#====================================================================
# 'simData' object constructor (initialization)
#====================================================================

simData <- function(desc, species, nsamples, initial.fit, setup.args, burnin,                     
                    nmcmc, l.int, Int, Slp, landings,  nHaul, nseas){
                    
    if (missing(species)|all(is.na(species))) stop("Missing 'species' parameter!!")
    if (missing(initial.fit))                 stop("Initital bayesian fit 'initial.fit' is missing!!")
    if (missing(l.int))                       stop("l.int is missing!!")
    if (missing(Int))                         stop("Int is missing!!")
    if (missing(Slp))                         stop("Slp is missing!!")
    if (missing(nHaul))                       stop("nHaul is missing!!")
    if (missing(nseas))                       stop("nseas is missing!!")
    if (missing(setup.args))                  stop("setup.args is missing!!")
    
    if(missing(nsamples)) nsamples <- 1
    if(missing(nmcmc))    nmcmc    <- 200
    if(missing(burnin))   burnin   <- nmcmc
    
    samples <- lapply(1:nsamples, function(x) costData())

    new("simData", desc = desc, species = species, samples = samples, initial.fit = initial.fit, 
                   setup.args = setup.args, burnin = burnin, nmcmc = nmcmc, l.int = l.int, 
                   Int = Int, Slp = Slp, landings = landings,  nHaul = nHaul, nseas = nseas)
 	  }

# simDataVal--------------------------------------------------------------------
setGeneric("simDataVal", function(obj,...){
	standardGeneric("simDataVal")
})
setMethod("simDataVal", signature("simData"), function(obj, desc,...){ 
        res <- new('simDataVal')
        slots <- slotNames(res)
        for(sl in slots) slot(res, sl) <- slot(obj, sl)
        res@desc <- ifelse(missing(desc), obj@desc,desc) 
        for(i in 1:length(obj@samples)) res@samples[[i]] <- costDataVal(obj@samples[[i]])
        return(res)})

# simDataCons-------------------------------------------------------------------
setGeneric("simDataCons", function(obj,objStrat,...){
	standardGeneric("simDataCons")
})

setMethod("simDataCons", signature("simDataVal", "strIni"), function(obj,
                                                                  objStrat,
                                                                  desc="Consolidated data",  
                                                                  ...){
        res <- new('simDataCons')
        slots <- slotNames(res)
        for(sl in slots) slot(res, sl) <- slot(obj, sl)
        res@desc <- ifelse(missing(desc), obj@desc,desc) 
        for(i in 1:length(obj@samples)) res@samples[[i]] <- costDataCons(obj = obj@samples[[i]], objStrat = objStrat, desc = desc,...)   
        
        return(res)                                                                
})    
                                  