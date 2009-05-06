#-------------------------------------------------------------------------------
# PerformStats class and creators.
# Dorleta Garcia :: Azti-Tecnalia
# 29/04/2009 15:50:50
#-------------------------------------------------------------------------------


setClass("PerformStats",
	representation(
    desc="character",                        # descriptor
    species="character",                     # recall of SL$spp (+ SL$sex)
	catchCat="character",                    # recall of the catch category (discards/landings)
	param="character",                       # recall of the parameter estimated (N, W, maturity, sex-ratio,...)
	strataDesc="strIni",                     # time, space and technical stratification considered
	methodDesc="character",                  # recall of the method (analytical, bootstrap, bayesian)
	nSamples = "numeric",                    # number of samples
    ageTrue = "data.frame",                  # True age structure (param-at-length)
    ageEst = "data.frame",                   # mean estimates of the age structure (param-at-length)
	ageAcc = "list",                         # a list with Accuracy statistics at age.
	ageBias = "list",                        # a list with Bias statistics at age.
	agePrec = "list",                        # a list with Precision statistics at age.
    lenTrue = "data.frame",               # True age structure (param-at-length)
    lenEst = "data.frame",                # mean estimates of the length structure (param-at-length)
	lenAcc = "list",                      # a list with Accuracy statistics at length.
	lenBias = "list",                     # a list with Bias statistics at length.
	lenPrec = "list",                     # a list with Precision statistics at length.
    totalNTrue = "data.frame",               # True age structure (param-at-length)
    totalNEst = "data.frame",                # mean estimates of the length structure (param-at-length)
	totalNAcc = "list",                      # a list with Accuracy statistics at length.
	totalNBias = "list",                     # a list with Bias statistics at length.
	totalNPrec = "list",                     # a list with Precision statistics at length.
	totalWTrue = "data.frame",               # True age structure (param-at-length)
    totalWEst = "data.frame",                # mean estimates of the length structure (param-at-length)
	totalWAcc = "list",                      # a list with Accuracy statistics at length.
	totalWBias = "list",                     # a list with Bias statistics at length.
	totalWPrec = "list"                     # a list with Precision statistics at length.
	),
	prototype(
        desc        = "PerformStats Object",
		species     = as.character(NA),
		catchCat    = as.character(NA),
		strataDesc  = strIni(),
		methodDesc  = as.character(NA),
		nSamples    = as.integer(NA),
# AGE -----------------------------------------------------------------------
        ageTrue     = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                age       = as.character(NA),
                                value     = as.numeric(NA)),
        ageEst      = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                age       = as.character(NA),
                                value     = as.numeric(NA)),
        ageAcc    = list(mse  = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           age       = as.character(NA),
                                           value     = as.numeric(NA)),
                         mae = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           age       = as.character(NA),
                                           value     = as.numeric(NA)),
                         smse = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           age       = as.character(NA),
                                           value     = as.numeric(NA)),
                         smae = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           age       = as.character(NA),
                                           value     = as.numeric(NA))),
        ageBias   = list(me  = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           age       = as.character(NA),
                                           value     = as.numeric(NA)),
                         sme = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           age       = as.character(NA),
                                           value     = as.numeric(NA)),
                         par = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           age       = as.character(NA),
                                           value     = as.numeric(NA)),
                         poe = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           age       = as.character(NA),
                                           value     = as.numeric(NA))),
        agePrec   = list(var  = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           age       = as.character(NA),
                                           value     = as.numeric(NA)),
                         cv = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           age       = as.character(NA),
                                           value     = as.numeric(NA))),
# LENGTH -----------------------------------------------------------------------
        lenTrue     = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                length       = as.character(NA),
                                value     = as.numeric(NA)),
        lenEst      = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                length       = as.character(NA),
                                value     = as.numeric(NA)),
        lenAcc    = list(mse  = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           length       = as.character(NA),
                                           value     = as.numeric(NA)),
                        mae = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           length       = as.character(NA),
                                           value     = as.numeric(NA)),
                         smse = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           length       = as.character(NA),
                                           value     = as.numeric(NA)),
                         smae = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           length       = as.character(NA),
                                           value     = as.numeric(NA))),
        lenBias   = list(me  = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           length       = as.character(NA),
                                           value     = as.numeric(NA)),
                         sme = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           length       = as.character(NA),
                                           value     = as.numeric(NA)),
                         par = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           length       = as.character(NA),
                                           value     = as.numeric(NA)),
                         poe = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           length       = as.character(NA),
                                           value     = as.numeric(NA))),
        lenPrec   = list(var  = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           length       = as.character(NA),
                                           value     = as.numeric(NA)),
                         cv = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           length       = as.character(NA),
                                           value     = as.numeric(NA))),
# total N -----------------------------------------------------------------------
        totalNTrue     = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                value     = as.numeric(NA)),
        totalNEst      = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                value     = as.numeric(NA)),
        totalNAcc    = list(mse  = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         mae = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         smse = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         smae = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA))),
        totalNBias   = list(me  = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         sme = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         par = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         poe = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA))),
        totalNPrec   = list(var  = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         cv = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA))),
# total W -----------------------------------------------------------------------
        totalWTrue     = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                value     = as.numeric(NA)),
        totalWEst      = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                value     = as.numeric(NA)),
        totalWAcc    = list(mse  = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         mae = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         smse = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         smae = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA))),
        totalWBias   = list(me  = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         sme = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         par = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         poe = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA))),
        totalWPrec   = list(var  = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         cv = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)))
                                           
    )
)

# Generic
setGeneric("PerformStats", function(estSimObj,
                                    trueDataObj,
                                 ...){
	standardGeneric("PerformStats")}
)


# 'missing' - 'missing'
setMethod("PerformStats", signature(estSimObj = "missing", trueDataObj ="missing"), function(estSimObj, trueDataObj,
                                                                                                             desc,...){
    res <- new('PerformStats', desc = desc)
    return(res)
})


# dbeOutputSim - trueData
setMethod("PerformStats", signature(estSimObj = "dbeOutputSim", trueDataObj ="trueData"), function(estSimObj, trueDataObj,
                                                                                                            desc, ...){
                                                                                                            
    if(!identical(estSimObj@strataDesc, trueDataObj@strataDesc)) 
        stop('The stratification in dbeOutputSim and trueData objects  must be the same!')
    
    desc <- ifelse(missing(desc), estSimObj@desc, desc)
 
    if(estSimObj@catchCat == 'LAN'){
        nSamples <- length(unique(estSimObj@lenStruc[['estim']]$sample))
        res      <- PS.lan(estSimObj, trueDataObj, desc = desc, nSamples = nSamples)
    }
    else{
        print('Function not available yet')
        res <- NULL
    }
    return(res)
  }
)


