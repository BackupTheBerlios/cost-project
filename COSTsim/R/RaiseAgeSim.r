#-------------------------------------------------------------------------------
# RaiseAgeSim and RaiseAgeSimBoot methods.
# Dorleta Garcia :: Azti-Tecnalia
# 29/04/2009 15:50:50
#-------------------------------------------------------------------------------


# RaiseAgeSim -----------------------------------------------------------------

setGeneric("RaiseAgeSim", function(dbeOutputSim,
                                 simObj,
                                 type="p",
                                 sex=as.character(NA),
                                 ...){
	standardGeneric("RaiseAgeSim")}
)



setMethod("RaiseAgeSim", signature(dbeOutputSim="dbeOutputSim",simObj="simDataCons"), function(dbeOutputSim,
                                                                                       simObj,
                                                                                       type="p",
                                                                                       sex=as.character(NA),
                                                                                       ...){
   nsamples <- length(simObj@samples)
   
   
   res1    <- dbeSimNULL(dbeOutputSim)
   res     <- dbeOutputSim
   
   res@ageStruc$estim <- res1@ageStruc$estim
   res@ageStruc$rep   <- res1@ageStruc$rep
   res@ageVar         <- res1@ageVar
   res@nSamp[['age']] <- res1@nSamp[['age']]
   res@nMeas[['age']] <- res1@nMeas[['age']]
   
  remove('res1')  
   
    for(i in 1:nsamples){
        dbeOut  <- dbeSim2dbe(dbeOutputSim, samples = i)
        xx <- RaiseAge(dbeOutput = dbeOut, csObject = simObj@samples[[i]]@cs,
                        type = type,  sex = sex, ...)

        slot(res,'nSamp')[['age']] <- rbind(slot(res,'nSamp')[['age']], cbind(sample = rep(i,dim(slot(xx,'nSamp')[['age']])[1]),slot(xx,'nSamp')[['age']]))
        slot(res,'nMeas')[['age']] <- rbind(slot(res,'nMeas')[['age']], cbind(sample = rep(i,dim(slot(xx,'nMeas')[['age']])[1]),slot(xx,'nMeas')[['age']]))

        d  <- ifelse(is.null(dim(slot(xx,'ageStruc')[['estim']])[1]), 1, dim(slot(xx,'ageStruc')[['estim']])[1])
        slot(res,'ageStruc')[['estim']] <- rbind(slot(res,'ageStruc')[['estim']], cbind(sample = rep(i,d),slot(xx,'ageStruc')[['estim']]))
        slot(res,'ageVar') <- rbind(slot(res,'ageVar'), cbind(sample = rep(i,d),slot(xx,'ageVar')))
    }

    return(res)
})



# RaiseAgeSimBoot --------------------------------------------------------------


setGeneric("RaiseAgeBootSim", function(dbeOutputSim,
                                 simObj,
                                 type="p",
                                 sex=as.character(NA),
                                 bootMethod = "samples",
                                 ...){
	standardGeneric("RaiseAgeBootSim")}
)


setMethod("RaiseAgeBootSim", signature(dbeOutputSim="dbeOutputSim",simObj="simDataCons"), function(dbeOutputSim,
                                                                                       simObj,
                                                                                       type="p",
                                                                                       sex=as.character(NA),
                                                                                       bootMethod = "samples",
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
        xx <- RaiseAgeBoot(dbeOutput = dbeOut, csObject = simObj@samples[[i]]@cs, type = 'p', 
                        sex=sex, bootMethod=bootMethod,...)
   
        slot(res,'nSamp')[['age']] <- rbind(slot(res,'nSamp')[['age']], cbind(sample = rep(i,dim(slot(xx,'nSamp')[['age']])[1]),slot(xx,'nSamp')[['age']]))
        slot(res,'nMeas')[['age']] <- rbind(slot(res,'nMeas')[['age']], cbind(sample = rep(i,dim(slot(xx,'nMeas')[['age']])[1]),slot(xx,'nMeas')[['age']]))

        d  <- ifelse(is.null(dim(slot(xx,'ageStruc')[['estim']])[1]), 1, dim(slot(xx,'ageStruc')[['estim']])[1])
        slot(res,'ageStruc')[['estim']] <- rbind(slot(res,'ageStruc')[['estim']], cbind(sample = rep(i,d),slot(xx,'ageStruc')[['estim']]))
        slot(res,'ageStruc')[['rep']]   <- rbind(slot(res,'ageStruc')[['rep']], cbind(sample = rep(i,dim(slot(xx,'ageStruc')[['rep']])[1]),slot(xx,'ageStruc')[['rep']]))
        slot(res,'ageVar')              <- rbind(slot(res,'ageVar'), cbind(sample = rep(i,d),slot(xx,'ageVar')))
    
    }
return(res)

})



