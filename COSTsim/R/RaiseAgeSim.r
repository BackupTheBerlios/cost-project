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
   
   
   res1    <- dbeSimNULL(dbeOutputSim)
   res     <- dbeOutputSim
   
   res@ageStruc$estim <- res1@ageStruc$estim
   res@ageStruc$rep   <- res1@ageStruc$rep
   res@ageVar         <- res1@ageVar
   
  remove('res1')  
   
    for(i in 1:nsamples){
        dbeOut  <- dbeSim2dbe(dbeOutputSim, samples = i)
        xx <- RaiseAge(dbeOutput = dbeOut, csObject = simObj@samples[[i]]@cs,
                        type = type,  sex = sex, ...)

        d  <- ifelse(is.null(dim(slot(xx,'ageStruc')[['estim']])[1]), 1, dim(slot(xx,'ageStruc')[['estim']])[1])
        slot(res,'ageStruc')[['estim']] <- rbind(slot(res,'ageStruc')[['estim']], cbind(sample = rep(i,d),slot(xx,'ageStruc')[['estim']]))
        slot(res,'ageVar') <- rbind(slot(res,'ageVar'), cbind(sample = rep(i,d),slot(xx,'ageVar')))
    
    
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
   
   dbeOut  <- dbeSim2dbe(dbeOutputSim)
   res1    <- dbeSimNULL(dbeOutputSim)
   res     <- dbeOutputSim
   
   res@ageStruc$estim <- res1@ageStruc$estim
   res@ageStruc$rep <- res1@ageStruc$rep
   res@ageVar         <- res1@ageVar
   
    remove('res1')
    for(i in 1:nsamples){
        xx <- RaiseAgeBoot(dbeOutput = dbeOut, csObject = simObj@samples[[i]]@cs, type = 'fixed', 
                        sex=sex, bootMethod=bootMethod,...)

        d  <- ifelse(is.null(dim(slot(xx,'ageStruc')[['estim']])[1]), 1, dim(slot(xx,'ageStruc')[['estim']])[1])
        slot(res,'ageStruc')[['estim']] <- rbind(slot(res,'ageStruc')[['estim']], cbind(sample = rep(i,d),slot(xx,'ageStruc')[['estim']]))
        slot(res,'ageStruc')[['rep']] <- rbind(slot(res,'ageStruc')[['rep']], cbind(sample = rep(i,d),slot(xx,'ageStruc')[['rep']]))
        slot(res,'ageVar') <- rbind(slot(res,'ageVar'), cbind(sample = rep(i,d),slot(xx,'ageVar')))
    
    }
return(res)

})



