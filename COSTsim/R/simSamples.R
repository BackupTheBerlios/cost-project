#-------------------------------------------------------------------------------
# simSamples method.
# Dorleta Garcia :: Azti-Tecnalia
# 29/04/2009 15:50:50
#-------------------------------------------------------------------------------

setGeneric("simSamples", function(obj,
                                 ndisc,
                                 ntrip,
                                 fit){
	standardGeneric("simSamples")}
)


setMethod("simSamples", signature(obj = "simDataCons"), function(obj, 
                                                                 ndisc = simDataObj@setup.args$ndisc, 
                                                                 ntrip = simDataObj@setup.args$ntrip, 
                                                                 fit = FALSE){
    simDataObj <- obj
    params     <- simDataObj@initial.fit
    class(params) <- 'fit.caa'
    nmcmc      <- simDataObj@nmcmc
    l.int      <- simDataObj@l.int
    Int        <- simDataObj@Int
    Slp        <- simDataObj@Slp
    landings   <- simDataObj@landings
    nHaul      <- simDataObj@nHaul
    nseas      <- simDataObj@nseas
    nsamples   <- length(simDataObj@samples)
    setup.data       <- simDataObj@setup.args
    setup.data$ndisc <- ndisc
    setup.data$ntrip <- ntrip
    setup.data       <- setup(setup.data)
    
    burnin     <- simDataObj@burnin
    
    fit.by <- vector('list', nsamples)

    
    for(i in 1:nsamples){
    cat('Sample: ', i,'\n')
          sim.out <- cost.simloop(params,setup.data, burnin,nmcmc,l.int,Int,Slp,
                        landings,nHaul,nseas, fit = fit)
          slot(simDataObj, 'samples')[[i]] <- convert2cost(sim.out)
          fit.by[[i]] <- sim.out$mbe.fit
    }

    if(fit == TRUE)
        return(list(simData = simDataObj, mbe.fit = fit.by))
    else
        return(simDataObj)
}
)