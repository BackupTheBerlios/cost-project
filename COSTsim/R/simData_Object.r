#-------------------------------------------------------------------------------
# simData, simDataVal, simDataCons classes and creators
# Dorleta Garcia :: Azti-Tecnalia
# 29/04/2009 15:50:50
#-------------------------------------------------------------------------------

#library(COSTcore)
#
#====================================================================
# dbeObj : outcome object for all COSTdbe package methods ouputs
#====================================================================


setClass("simData",
	representation(
        desc        = "character",                     # descriptor
        species     = "character",                     # recall of SL$spp (+ SL$sex)
        samples     = "list",                          # list with simulates 'csData' objects.
        cl          = "clData",
        ce          = "ceData",
        initial.fit = "list",
		setup.args  = "list"                          # list with setup parameters
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
        samples     = "list",                          # list with simulates 'csData' objects.
        cl          = "clDataVal",
        ce          = "ceDataVal",
        initial.fit = "list",
		setup.args  = "list"                           # list with setup parameters
		))
		
setClass("simDataCons",
	representation(
        desc        = "character",                     # descriptor
        species     = "character",                     # recall of SL$spp (+ SL$sex)
        samples     = "list",                          # list with simulates 'csData' objects.
        cl          = "clDataCons",
        ce          = "ceDataCons",
        initial.fit = "list",
		setup.args  = "list"                                # list with setup parameters
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

# simData --------------------------------------------------------------------
setGeneric("simData", function(obj, ...){
	standardGeneric("simData")
})
setMethod("simData", signature(obj = 'missing'),   function(obj, desc, species, nsamples, 
                    initial.fit, setup.args, burnin,                     
                    nmcmc, l.int, Int, Slp, landings,  nHaul, nseas){ 
    res <- simData(desc, species, nsamples, initial.fit, setup.args, burnin,                     
                    nmcmc, l.int, Int, Slp, landings,  nHaul, nseas)
return(res)})


simData <- function(desc, species, nsamples, clObj, ceObj, initial.fit, setup.args){
                    
    if (missing(species)|all(is.na(species))) stop("Missing 'species' parameter!!")
    if (missing(initial.fit))                 stop("Initital bayesian fit 'initial.fit' is missing!!")
#    if (missing(l.int))                       stop("l.int is missing!!")
#    if (missing(Int))                         stop("Int is missing!!")
#    if (missing(Slp))                         stop("Slp is missing!!")
#    if (missing(nHaul))                       stop("nHaul is missing!!")
#    if (missing(nseas))                       stop("nseas is missing!!")
    if (missing(setup.args))                  stop("setup.args is missing!!")
    
    if(missing(nsamples)) nsamples <- 1

    samples <- lapply(1:nsamples, function(x) csData())

    new("simData", desc = desc, species = species, samples = samples, initial.fit = initial.fit, 
                   setup.args = setup.args, cl = clObj, ce = ceObj)
 	  }

# simDataVal--------------------------------------------------------------------
setGeneric("simDataVal", function(obj,...){
	standardGeneric("simDataVal")
})
setMethod("simDataVal", signature("simData"), function(obj, desc,...){ 
        res <- new('simDataVal')
        slots <- c("species", "initial.fit", "setup.args")
        for(sl in slots) slot(res, sl) <- slot(obj, sl)
        res@desc <- ifelse(missing(desc), obj@desc,desc) 
        res@cl  <- clDataVal(obj@cl)
        res@ce  <- ceDataVal(obj@ce)
        for(i in 1:length(obj@samples)) res@samples[[i]] <- csDataVal(obj@samples[[i]])
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
        slots <- c("species", "initial.fit", "setup.args")
        for(sl in slots) slot(res, sl) <- slot(obj, sl)
        res@desc <- ifelse(missing(desc), obj@desc,desc) 
        res@cl <- clDataCons(obj = obj@cl, objStrat = objStrat, ...) 
        res@ce <- ceDataCons(obj = obj@ce, objStrat = objStrat, ...)
        for(i in 1:length(obj@samples)) res@samples[[i]] <- csDataCons(obj = obj@samples[[i]], objStrat = objStrat, desc = desc,...)   
        
        return(res)                                                                
})    
                                  