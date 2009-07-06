#-------------------------------------------------------------------------------
# totVolumeSim method.
# Dorleta Garcia :: Azti-Tecnalia
# 29/04/2009 15:50:50
#-------------------------------------------------------------------------------

setGeneric("totVolumeSim", function(dbeOutputSim, simObj, ...){
	standardGeneric("totVolumeSim")}
)


setMethod("totVolumeSim", signature(dbeOutput="dbeOutputSim",simObj = "simDataCons"),
    function(dbeOutputSim, simObj, type="trip", val="weight", sampPar = TRUE, ...){
    
   nsamples <- length(simObj@samples)

   bs  <- names(which(!(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x))) %in% c('list', 'data.frame'))))
   df  <- names(which(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x)))== 'data.frame'))
   lst <- names(which(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x)))== 'list'))
   lst.df <- sapply(lst, function(x) lapply(slot(dbeObject(species = '1'),x), class))


   dbeOut <- dbeSim2dbe(dbeOutputSim,1)
   res    <- dbeSimNULL(dbeOutputSim,1)

   if(type == 'landings'){
        for(i in 1:nsamples){
            xx <- totVolume(dbeOutput = dbeOut, csObject = simObj@samples[[i]],
                                                ceObject = simObj@ce,
                                                clObject = simObj@cl,
                        type = type,  val = val, sampPar = sampPar,...)

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
    }
    
    else{ # type != landings
        for(i in 1:nsamples){
            xx <- totVolume(dbeOutput = dbeOut, csObject = simObj@samples[[i]],
                                                ceObject = simObj@ce,
                        type = type,  val = val, sampPar = sampPar,...)
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
    }
    
    return(res)

})


tvs <- function(dbeOutputSim, simObj, type="trip", val="weight", sampPar = TRUE, ...){
    
   nsamples <- length(simObj@samples)

   bs  <- names(which(!(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x))) %in% c('list', 'data.frame'))))
   df  <- names(which(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x)))== 'data.frame'))
   lst <- names(which(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x)))== 'list'))
   lst.df <- sapply(lst, function(x) lapply(slot(dbeObject(species = '1'),x), class))


   dbeOut <- dbeSim2dbe(dbeOutputSim,1)
   res    <- dbeSimNULL(dbeOutputSim)

   if(type == 'landings'){
        for(i in 1:nsamples){
            xx <- totVolume(dbeOutput = dbeOut, csObject = simObj@samples[[i]],
                                                ceObject = simObj@ce,
                                                clObject = simObj@cl,
                        type = type,  val = val, sampPar = sampPar,...)

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
    }
    
    else{ # type != landings
        for(i in 1:nsamples){
            xx <- totVolume(dbeOutput = dbeOut, csObject = simObj@samples[[i]],
                                                ceObject = simObj@ce,
                        type = type,  val = val, sampPar = sampPar,...)
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
    }
    
    return(res)

}



