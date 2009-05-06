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
                                # if TRUE, updated 'dbeOutput' object is returned ; if FALSE, only resulting dataframe is returned
                               ...){
standardGeneric("dbeCalcSim")
})

 
setMethod("dbeCalcSim",signature(object="dbeOutputSim"),function(object,              #'dbeOutput' object 
                                                           type="CI",           # "CI" for confidence interval, or "CV" for coefficient of variation
                                                           vrbl="l",            # specifies data on which calculation is applied : "l" for length structure, "a" for age structure, "n" for total number estimates, "w" for total weight estimates
                                                           probs=c(0.025,0.975),# used only if type="CI", defines bounds
                                                           replicates=FALSE,    # if TRUE, calculation is made from $rep elements ; if FALSE, $estim and @...Var are used 
                                                           ...){                                

    outpuT <- NULL
    #according to input parameters, one of the previous functions is to be used
    lType <- tolower(type)
    if (!lType%in%c("ci","cv")) stop("wrong 'type' parameter!!") 
    if (!is.logical(replicates)) stop("wrong 'replicates' parameter!!") 
    
    sl <- switch(vrbl,
                l = 'lenNum',
                a = 'ageNum',
                n= 'totalNnum',
                w= 'totalWnum')
    dat <- switch(vrbl,
                l = 'lenStruc',
                a = 'ageStruc',
                n= 'totalN',
                w= 'totalW')
                
    df <- switch(type,
                CI = 'ci',
                CV = c('cv', 'DCRcvIndicator'))

    res <- object
    slot(res, sl)[[df[1]]] <- NULL
    
    nsamples <- length(unique(slot(object, dat)[['estim']]$sample))
    
    
    if(df[1] == 'cv')  slot(res, sl)[[df[2]]] <- NULL

    for(i in 1:nsamples){

        dbeOut <- dbeSim2dbe(object, samples = i)
        xx <- dbeCalc(object = dbeOut, type = type, vrbl = vrbl, probs = probs, replicates = replicates, update = FALSE, ...)
        
        d  <- ifelse(class(xx) == 'list', dim(xx[[1]])[1], dim(xx)[1])
          
        if(df[1] == 'ci'){
            slot(res,sl)[[df[1]]] <- rbind(slot(res,sl)[[df[1]]], cbind(sample = rep(i,d), xx))
        }
        else{ # CV
            slot(res,sl)[[df[1]]] <- rbind(slot(res,sl)[[df[1]]], cbind(sample = rep(i,d), xx[['DF']]))
            slot(res,sl)[[df[2]]] <- as.data.frame(rbind(slot(res,sl)[[df[2]]], cbind(sample = i, value = xx[['dcrInd']])))
        }
    }
        
    
 return(res)

})





                                                                   