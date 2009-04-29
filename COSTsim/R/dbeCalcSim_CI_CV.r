#-------------------------------------------------------------------------------
# dbeCalcSim method.
# Dorleta Garcia :: Azti-Tecnalia
# 29/04/2009 15:50:50
#-------------------------------------------------------------------------------


        
setGeneric("dbeCalcSim", function(object,              # 'dbeOutput' object
                               type="CI",           # "CI" for confidence interval, or "CV" for coefficient of variation
                               vrbl="l",            # specifies data on which calculation is applied : "l" for length structure, "a" for age structure, "n" for total number estimates, "w" for total weight estimates
                               probs=c(0.025,0.975),# used only if type="CI", defines bounds
                               replicates=FALSE,    # if TRUE, calculation is made from $rep elements ; if FALSE, $estim and @...Var are used 
                               update=FALSE,        # if TRUE, updated 'dbeOutput' object is returned ; if FALSE, only resulting dataframe is returned
                               ...){
standardGeneric("dbeCalcSim")
})



 
setMethod("dbeCalcSim",signature(object="dbeOutputSim"),function(object,              #'dbeOutput' object 
                                                           type="CI",           # "CI" for confidence interval, or "CV" for coefficient of variation
                                                           vrbl="l",            # specifies data on which calculation is applied : "l" for length structure, "a" for age structure, "n" for total number estimates, "w" for total weight estimates
                                                           probs=c(0.025,0.975),# used only if type="CI", defines bounds
                                                           replicates=FALSE,    # if TRUE, calculation is made from $rep elements ; if FALSE, $estim and @...Var are used 
                                                           update=FALSE,        # if TRUE, updated 'dbeOutput' object is returned ; if FALSE, only resulting dataframe is returned
                                                           ...){                                

    outpuT <- NULL
    #according to input parameters, one of the previous functions is to be used
    lType <- tolower(type)
    if (!lType%in%c("ci","cv")) stop("wrong 'type' parameter!!") 
    if (!is.logical(replicates)) stop("wrong 'replicates' parameter!!") 
    if (!is.logical(update)) stop("wrong 'update' parameter!!")

    nsamples <- length(unique(object@nSamp$len$sample))
    
    bs  <- names(which(!(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x))) %in% c('list', 'data.frame'))))
    df  <- names(which(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x)))== 'data.frame'))
    lst <- names(which(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x)))== 'list'))
    lst.df <- sapply(lst, function(x) lapply(slot(dbeObject(species = '1'),x), class))


   res    <- dbeSimNULL(object)

    for(i in 1:nsamples){
    
    
        dbeOut <- dbeSim2dbe(object, samples = i)
        xx <- dbeCalc(object = dbeOut, type = type, vrbl = vrbl, probs = probs, replicates = replicates, update = update, ...)
        
        for(s in df){
            d  <- dim(slot(xx,s))[1]
            slot(res,s) <- rbind(slot(res,s), cbind(sample = rep(i,d),slot(xx,s)))
        }
        for(s in lst){
            for(sl in names(lst.df[[s]])){
                d  <- ifelse(is.null(dim(slot(xx,s)[[sl]])[1]), 1, dim(slot(xx,s)[[sl]])[1])
                slot(res,s)[[sl]] <- rbind(slot(res,s)[[sl]], cbind(sample = rep(i,d),slot(xx,s)[[sl]]))
            }
        }
    }
 return(res)

})



 
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

## 'stratAgreg' method for aggegating 'dbeOutput' tables (estim & var) : other output tables are empty
## WARNING : summing the variances requires strong probability assumptions
 

setGeneric("stratAggreg", function(object,                  # 'dbeOutput' object
                                  timeStrata=TRUE,         # if TRUE, aggregation is made over time strata
                                  spaceStrata=TRUE,        # if TRUE, aggregation is made over space strata
                                  techStrata=FALSE,        # if TRUE, aggregation is made over technical strata
                                  ...){
standardGeneric("stratAggreg")
})
setMethod("stratAggreg", signature(object="dbeOutputSim"),function(object,                  # 'dbeOutput' object
                                                               timeStrata=TRUE,         # if TRUE, aggregation is made over time strata
                                                               spaceStrata=TRUE,        # if TRUE, aggregation is made over space strata
                                                               techStrata=FALSE,        # if TRUE, aggregation is made over technical strata
                                                               ...){

    #'strataDesc' field is updated according to input parameters                                              #
    if (timeStrata) {object@strataDesc@timeStrata <- NA ; object@strataDesc@tpRec <- list(NA)}                # ADDED MM : 22/04/2009 
    if (spaceStrata) {object@strataDesc@spaceStrata <- NA ; object@strataDesc@spRec <- list(NA)}              #
    if (techStrata) {object@strataDesc@techStrata <- NA ; object@strataDesc@tcRec <- list(NA)}                #

    #subfunction applied to each table
    agg <- function(tab,nSampAge=FALSE) {

        if (all(is.na(tab))) {
            return(tab)
        } else {
            if (timeStrata) tab$time <- "all"
            if (spaceStrata) tab$space <- "all"
            if (!nSampAge & techStrata) tab$technical <- "all"
                newTab <- aggregate(tab$value,as.list(tab[,(ncol(tab)-1):1]),sum,na.rm=TRUE)
                    names(newTab)[ncol(newTab)] <- "value"
        return(newTab[,names(tab)])}}
        

    object@nSamp$len <- agg(object@nSamp$len)
object@nSamp$age <- agg(object@nSamp$age,nSampAge=TRUE)
object@nMeas$len <- agg(object@nMeas$len)
object@nMeas$age <- agg(object@nMeas$age,nSampAge=TRUE)

object@lenStruc$estim <- agg(object@lenStruc$estim)
object@lenStruc$rep <- new("dbeOutput")@lenStruc$rep
object@lenVar <- agg(object@lenVar)
object@lenNum$ci <- new("dbeOutput")@lenNum$ci
object@lenNum$cv <- new("dbeOutput")@lenNum$cv
object@lenNum$DCRcvIndicator <- new("dbeOutput")@lenNum$DCRcvIndicator

object@ageStruc$estim <- agg(object@ageStruc$estim)
object@ageStruc$rep <- new("dbeOutput")@ageStruc$rep
object@ageVar <- agg(object@ageVar)
object@ageNum$ci <- new("dbeOutput")@ageNum$ci
object@ageNum$cv <- new("dbeOutput")@ageNum$cv
object@ageNum$DCRcvIndicator <- new("dbeOutput")@ageNum$DCRcvIndicator

object@totalN$estim <- agg(object@totalN$estim)
object@totalN$rep <- new("dbeOutput")@totalN$rep
object@totalNvar <- agg(object@totalNvar)
object@totalNnum$ci <- new("dbeOutput")@totalNnum$ci
object@totalNnum$cv <- new("dbeOutput")@totalNnum$cv
object@totalNnum$DCRcvIndicator <- new("dbeOutput")@totalNnum$DCRcvIndicator

object@totalW$estim <- agg(object@totalW$estim)
object@totalW$rep <- new("dbeOutput")@totalW$rep
object@totalWvar <- agg(object@totalWvar)
object@totalWnum$ci <- new("dbeOutput")@totalWnum$ci
object@totalWnum$cv <- new("dbeOutput")@totalWnum$cv
object@totalWnum$DCRcvIndicator <- new("dbeOutput")@totalWnum$DCRcvIndicator

return(object)

})


                                                                   