
setGeneric("totVolumeSim", function(dbeOutputSim, simObj, ...){
	standardGeneric("totVolumeSim")}
)


setMethod("totVolumeSim", signature(dbeOutput="dbeOutputSim",simObj = "simDataCons"),
    function(dbeOutputSim, simObj, type="trip", val="weight", sampPar = TRUE, rtl = FALSE, ...){
    
   nsamples <- length(simObj@samples)

   bs  <- names(which(!(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x))) %in% c('list', 'data.frame'))))
   df  <- names(which(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x)))== 'data.frame'))
   lst <- names(which(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x)))== 'list'))
   lst.df <- sapply(lst, function(x) lapply(slot(dbeObject(species = '1'),x), class))


   dbeOut <- dbeSim2dbe(dbeOutputSim)
   res    <- dbeSimNULL(dbeOutputSim)

   if(rtl == TRUE){
        for(i in 1:nsamples){
            xx <- totVolume(dbeOutput = dbeOut, csObject = simObj@samples[[i]]@cs,
                                                ceObject = simObj@samples[[i]]@ce,
                                                clObject = simObj@samples[[i]]@cl,
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
    
    else{ # rtl = FALSE
        for(i in 1:nsamples){
            xx <- totVolume(dbeOutput = dbeOut, csObject = simObj@samples[[i]]@cs,
                                                ceObject = simObj@samples[[i]]@ce,
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



