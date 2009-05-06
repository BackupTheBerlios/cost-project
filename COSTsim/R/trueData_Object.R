#-------------------------------------------------------------------------------
# trueData class and creators.
# Dorleta Garcia :: Azti-Tecnalia
# 29/04/2009 15:50:50
#-------------------------------------------------------------------------------

setClass("trueData",
	representation(
    desc       = "character",                 # descriptor
    species    = "character",                 # recall of SL$spp (+ SL$sex)
   	strataDesc = "strIni",                    # time, space and technical stratification considered
    lal        = "data.frame",                # True  landings at age
    laa        = "data.frame",                # True  landings at length
    ltw        = "data.frame",                # True  landings at age
    ltn        = "data.frame",                # True  landings at length
    dal        = "data.frame",                # True discards at length
    daa        = "data.frame",                 # True discards total weight
    dtw        = "data.frame",                # True  landings at age
    dtn        = "data.frame"                # True  landings at length
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
        ltw     = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                value     = as.numeric(NA)),
        ltn     = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                value     = as.numeric(NA)),
        dal     = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                length = as.character(NA),
                                value     = as.numeric(NA)),
        daa      = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                age = as.character(NA),
                                value     = as.numeric(NA)),
        dtw      = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                value     = as.numeric(NA)),     
        dtn      = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                value     = as.numeric(NA))                        
                                    )
)

setGeneric("trueData", function(obj, desc, ...){
	standardGeneric("trueData")}
)


setMethod("trueData", signature(obj = "missing"), function(obj, desc, species, strataDesc = strIni(), dal, daa, dtw, dtn, laa, lal, ltw, ltn){
    res <- new('trueData', desc = desc, species = species, strataDesc = strataDesc)
    if(!missing(dal)) res@dal <- dal
    if(!missing(daa)) res@daa <- daa
    if(!missing(dtw)) res@dtw <- dtw
    if(!missing(dtn)) res@dtn <- dtn
    if(!missing(lal)) res@lal <- lal
    if(!missing(laa)) res@laa <- laa
    if(!missing(ltw)) res@ltw <- ltw
    if(!missing(ltn)) res@ltn <- ltn
    
    return(res)
})




