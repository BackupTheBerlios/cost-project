#-------------------------------------------------------------------------------
# trueData class and creators.
# Dorleta Garcia :: Azti-Tecnalia
# 29/04/2009 15:50:50
#-------------------------------------------------------------------------------

setClass("trueData",
	representation(
    desc       = "character",                 # descriptor
    species    = "character",                 # recall of SL$spp (+ SL$sex)
   	strataDesc = "strIni",                     # time, space and technical stratification considered
    lal        = "data.frame",                # True  landings at age
    laa        = "data.frame",                # True  landings at length
    dal        = "data.frame",                # True discards at length
    dtw        = "data.frame"                # True discards total weight
	),
	prototype(
        desc        = "PerformStats Object",
		species     = as.character(NA),
		strataDesc  = strIni(),
        lal     = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                length = as.character(NA),
                                value     = as.numeric(NA)),
        laa      = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                age = as.character(NA),
                                value     = as.numeric(NA)),
        dal     = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                length = as.character(NA),
                                value     = as.numeric(NA)),
        dtw      = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                value     = as.numeric(NA))                        
                                    )
)

setGeneric("trueData", function(obj, desc, ...){
	standardGeneric("trueData")}
)


setMethod("trueData", signature(obj = "missing"), function(obj, desc, species, strataDesc = strIni(), dal, dtw, laa, lal){
    res <- new('trueData', desc = desc, species = species, strataDesc = strataDesc)
    if(!missing(dal)) res@dal <- dal
    if(!missing(dtw)) res@dtw <- dtw
    if(!missing(laa)) res@dal <- laa
    if(!missing(lal)) res@dal <- lal
    
    return(res)
})




