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


setClass("mbeOutputSim",
	representation(
    LAN = "OutputSim",                        #descriptor
    DIS = "OutputSim"),                     #recall of SL$spp (+ SL$sex)
	
    prototype(
        LAN     = ObjectSim(species = 'zzz'),
		DIS     = ObjectSim(species = 'zzz'))
)



#====================================================================
# 'dbeOutput' object constructor (initialization)
#====================================================================

mbeObjectSim <- function(desc, species, strataDesc){

    if (missing(desc)) desc <- as.character(NA)
    if (missing(species)|all(is.na(species))) stop("Missing 'species' parameter!!")
    if (missing(strataDesc)) strataDesc <- strIni()

    res <- new("mbeOutputSim")
    res@LAN <- new("OutputSim",desc=desc,species=species,catchCat= 'LAN',strataDesc=strataDesc, methodDesc='bayesian')
    res@DIS <- new("OutputSim",desc=desc,species=species,catchCat= 'DIS',strataDesc=strataDesc, methodDesc='bayesian')
    return(res)
 	  }


