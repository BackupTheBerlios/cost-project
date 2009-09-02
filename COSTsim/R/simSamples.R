#-------------------------------------------------------------------------------
# simSamples method.
# Dorleta Garcia :: Azti-Tecnalia
# 29/04/2009 15:50:50
#-------------------------------------------------------------------------------

setGeneric("simSamples", function(obj,
                                ...){
	standardGeneric("simSamples")}
)


#setup.data<-new.setup(seas.list,area.list,gear.list,species,cov.list,ntrip=n.ml.trips+n.observer.trips,
#                      ndisc=n.observer.trips,div=c(4,9),ageMin,ageMax)




setMethod("simSamples", signature(obj = "simData"), function(obj, n.observer.trips, n.ml.trips,
                                                                 nlsamp.land, nasamp.land, nlsamp.disc, nasamp.disc){

    obj@setup.args$n.observer.trips  <- ifelse(missing(n.observer.trips), obj@setup.args$n.observer.trips, n.observer.trips)
    obj@setup.args$n.ml.trips        <- ifelse(missing(n.ml.trips), obj@setup.args$n.ml.trips, n.ml.trips)
    obj@setup.args$nlsamp.land       <- ifelse(missing(nlsamp.land), obj@setup.args$nlsamp.land, nlsamp.land)
    obj@setup.args$nasamp.land       <- ifelse(missing(nasamp.land), obj@setup.args$nasamp.land, nasamp.land)
    obj@setup.args$nlsamp.disc       <- ifelse(missing(nlsamp.disc), obj@setup.args$nlsamp.disc, nlsamp.disc)
    obj@setup.args$nasamp.disc       <- ifelse(missing(nasamp.disc), obj@setup.args$nasamp.disc, nasamp.disc)

   # setup.data       <- setup(setup.args)


    for(i in 1:length(obj@samples)){
        cat('Sample: ', i,'\n')
        setup.data<-new.setup(obj@initial.fit$fit, obj@setup.args$use.seasons, obj@setup.args$arealist, obj@setup.args$use.gears, obj@setup.args$species,
                        obj@setup.args$age.covariates,obj@setup.args$weight.covariates, nmland= obj@setup.args$n.ml.trips,
                      nobs=obj@setup.args$n.observer.trips,obj@setup.args$ageMin,obj@setup.args$ageMax)
        sim.out <- make.sim.data(obj@initial.fit$fit,setup.data,  obj@setup.args$nlsamp.land,  obj@setup.args$nasamp.land,  obj@setup.args$nlsamp.disc,
                     obj@setup.args$nasamp.disc,  obj@setup.args$length.list)

slot(obj, 'samples')[[i]] <- convert2cost(sim.out, species = obj@setup.args$species, gear.list = setup.data$gearlist, 
    area.list= setup.data$arealist, year = obj@setup.args$year)

    }
    return(obj)
}
)