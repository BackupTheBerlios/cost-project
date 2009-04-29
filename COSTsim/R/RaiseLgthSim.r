#-------------------------------------------------------------------------------
# RaiseLgthSim and RaiseLgthSimBoot methods.
# Dorleta Garcia :: Azti-Tecnalia
# 29/04/2009 15:50:50
#-------------------------------------------------------------------------------


# RaiseLgthSim -----------------------------------------------------------------

setGeneric("RaiseLgthSim", function(dbeOutputSim,
                                 simData,
                                 spp,
                                 taxon,
                                 sex=as.character(NA),
                                 ...){
	standardGeneric("RaiseLgthSim")}
)


setMethod("RaiseLgthSim", signature(dbeOutputSim = "dbeOutputSim", simData ="simDataCons"), function(dbeOutputSim, simData,
                                                                                                              spp,
                                                                                                              taxon,
                                                                                                              sex=as.character(NA),...){
   nsamples <- length(simData@samples)
   
   bs  <- names(which(!(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x))) %in% c('list', 'data.frame'))))
   df  <- names(which(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x)))== 'data.frame'))
   lst <- names(which(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x)))== 'list'))
   lst.df <- sapply(lst, function(x) lapply(slot(dbeObject(species = '1'),x), class))


   dbeOut <- dbeSim2dbe(dbeOutputSim)
   res    <- dbeSimNULL(dbeOutputSim)

    for(i in 1:nsamples){
        xx <- RaiseLgth(dbeOutput = dbeOut, csObject = simData@samples[[i]]@cs, clObject = simData@samples[[i]]@cl,
                        spp = spp, taxon = taxon, sex = sex,...)

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


# RaiseLgthSimBoot--------------------------------------------------------------

setGeneric("RaiseLgthBootSim", function(dbeOutputSim,
                                 simData,
                                 spp,
                                 taxon,
                                 sex=as.character(NA),
                                 B,
                                 ...){
	standardGeneric("RaiseLgthBootSim")}
)

setMethod("RaiseLgthBootSim", signature(dbeOutputSim="dbeOutputSim",simData="simDataCons"), function(dbeOutputSim,simData,
                                                                                                              spp,
                                                                                                              taxon,
                                                                                                              sex=as.character(NA),
                                                                                                              B,
                                                                                                              ...){
                                                                                                           
   nsamples <- length(simData@samples)
   
   bs  <- names(which(!(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x))) %in% c('list', 'data.frame'))))
   df  <- names(which(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x)))== 'data.frame'))
   lst <- names(which(sapply(slotNames(dbeObject(species = '1')), function(x) class(slot(dbeObject(species = '1'), x)))== 'list'))
   lst.df <- sapply(lst, function(x) lapply(slot(dbeObject(species = '1'),x), class))


   dbeOut <- dbeSim2dbe(dbeOutputSim)
   res    <- dbeSimNULL(dbeOutputSim)

    for(i in 1:nsamples){
        xx <- RaiseLgthBoot(dbeOutput = dbeOut, csObject = simData@samples[[i]]@cs, clObject = simData@samples[[i]]@cl,
                        spp = spp, taxon = taxon, sex = sex, B = B,...)

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
