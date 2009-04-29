#-------------------------------------------------------------------------------
# RaiseAgeSim and RaiseAgeSimBoot methods.
# Dorleta Garcia :: Azti-Tecnalia
# 29/04/2009 15:50:50
#-------------------------------------------------------------------------------


# RaiseAgeSim -----------------------------------------------------------------

setGeneric("RaiseAgeSim", function(dbeOutputSim,
                                 simObj,
                                 type="fixed",
                                 sex=as.character(NA),
                                 ...){
	standardGeneric("RaiseAgeSim")}
)



setMethod("RaiseAgeSim", signature(dbeOutputSim="dbeOutputSim",simObj="simDataCons"), function(dbeOutputSim,
                                                                                       simObj,
                                                                                       type="fixed",
                                                                                       sex=as.character(NA),
                                                                                       ...){
   nsamples <- length(simObj@samples)
   
   bs  <- names(which(!(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x))) %in% c('list', 'data.frame'))))
   df  <- names(which(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x)))== 'data.frame'))
   lst <- names(which(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x)))== 'list'))
   lst.df <- sapply(lst, function(x) lapply(slot(dbeObject(species = '1'),x), class))


   dbeOut <- dbeSim2dbe(dbeOutputSim)
   res    <- dbeSimNULL(dbeOutputSim)

    for(i in 1:nsamples){
        xx <- RaiseAge(dbeOutput = dbeOut, csObject = simObj@samples[[i]]@cs,
                        type = type,  sex = sex,...)

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



# RaiseAgeSimBoot --------------------------------------------------------------


setGeneric("RaiseAgeBootSim", function(dbeOutputSim,
                                 simObj,
                                 type="fixed",
                                 sex=as.character(NA),
                                 bootMethod = "samples",
                                 ...){
	standardGeneric("RaiseAgeBootSim")}
)


setMethod("RaiseAgeBootSim", signature(dbeOutputSim="dbeOutputSim",simObj="simDataCons"), function(dbeOutputSim,
                                                                                       simObj,
                                                                                       type="fixed",
                                                                                       sex=as.character(NA),
                                                                                       bootMethod = "samples",
                                                                                       ...){
                                                                                                           

    nsamples <- length(simObj@samples)
   
   bs  <- names(which(!(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x))) %in% c('list', 'data.frame'))))
   df  <- names(which(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x)))== 'data.frame'))
   lst <- names(which(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x)))== 'list'))
   lst.df <- sapply(lst, function(x) lapply(slot(dbeObject(species = '1'),x), class))


   dbeOut <- dbeSim2dbe(dbeOutputSim)
   res    <- dbeSimNULL(dbeOutputSim)

    for(i in 1:nsamples){
        xx <- RaiseAgeBoot(dbeOutput = dbeOut, csObject = simObj@samples[[i]]@cs, type = 'fixed', 
                        sex=sex, bootMethod=bootMethod,...)

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



