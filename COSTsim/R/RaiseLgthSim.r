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
   
   res1    <- dbeSimNULL(dbeOutputSim)
   res     <- dbeOutputSim
   
   res@lenStruc[['estim']]     <- res1@lenStruc$estim
   res@lenStruc[['rep']]       <- res1@lenStruc$rep
   res@lenVar                  <- res1@lenVar
   res@nSamp[['len']]          <- res1@nSamp[['len']]
   res@nMeas[['len']]          <- res1@nMeas[['len']]
     
    remove('res1')
    dbeOut  <- dbeSim2dbe(dbeOutputSim)
    for(i in 1:nsamples){
        
        xx <- RaiseLgth(dbeOutput = dbeOut, csObject = simData@samples[[i]], clObject = simData@cl,
                        spp = spp, taxon = taxon, sex = sex,...)
        
        slot(res,'nSamp')[['len']] <- rbind(slot(res,'nSamp')[['len']], cbind(sample = rep(i,dim(slot(xx,'nSamp')[['len']])[1]),slot(xx,'nSamp')[['len']]))
        slot(res,'nMeas')[['len']] <- rbind(slot(res,'nMeas')[['len']], cbind(sample = rep(i,dim(slot(xx,'nMeas')[['len']])[1]),slot(xx,'nMeas')[['len']]))

        d  <- ifelse(is.null(dim(slot(xx,'lenStruc')[['estim']])[1]), 1, dim(slot(xx,'lenStruc')[['estim']])[1])
        slot(res,'lenStruc')[['estim']] <- rbind(slot(res,'lenStruc')[['estim']], cbind(sample = rep(i,d),slot(xx,'lenStruc')[['estim']]))
        slot(res,'lenVar')              <- rbind(slot(res,'lenVar'), cbind(sample = rep(i,d),slot(xx,'lenVar')))

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
   
   res1    <- dbeSimNULL(dbeOutputSim)
   res     <- dbeOutputSim
   
   res@lenStruc$estim <- res1@lenStruc$estim
   res@lenStruc$rep   <- res1@lenStruc$rep
   res@lenVar         <- res1@lenVar
   res@nSamp[['len']]          <- res1@nSamp[['len']]
   res@nMeas[['len']]          <- res1@nMeas[['len']]
   
   dbeOut  <- dbeSim2dbe(dbeOutputSim)
    
  remove('res1')
    for(i in 1:nsamples){
        
        xx <- RaiseLgthBoot(dbeOutput = dbeOut, csObject = simData@samples[[i]], clObject = simData@cl,
                        spp = spp, taxon = taxon, sex = sex, B = B,...)
        
        slot(res,'nSamp')[['len']] <- rbind(slot(res,'nSamp')[['len']], cbind(sample = rep(i,dim(slot(xx,'nSamp')[['len']])[1]),slot(xx,'nSamp')[['len']]))
        slot(res,'nMeas')[['len']] <- rbind(slot(res,'nMeas')[['len']], cbind(sample = rep(i,dim(slot(xx,'nMeas')[['len']])[1]),slot(xx,'nMeas')[['len']]))

        d  <- ifelse(is.null(dim(slot(xx,'lenStruc')[['estim']])[1]), 1, dim(slot(xx,'lenStruc')[['estim']])[1])
        slot(res,'lenStruc')[['estim']] <- rbind(slot(res,'lenStruc')[['estim']], cbind(sample = rep(i,d),  slot(xx,'lenStruc')[['estim']]))
        slot(res,'lenStruc')[['rep']]   <- rbind(slot(res,'lenStruc')[['rep']],   cbind(sample = rep(i,dim(slot(xx,'lenStruc')[['rep']])[1]),slot(xx,'lenStruc')[['rep']]))
        slot(res,'lenVar')              <- rbind(slot(res,'lenVar'), cbind(sample = rep(i,d),slot(xx,'lenVar')))
    
    }
return(res)

})
