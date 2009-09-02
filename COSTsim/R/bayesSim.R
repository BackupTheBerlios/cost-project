# By now it only extracts results related to landings but it would be quite 
# straightforward (copy-paste and  change lan by dis) to put the discards.

bayesSim <- function(sim.dt, setup.args, nmcmc, thin, burnin, rep = TRUE){

        species <- sim.dt@species

        res <- mbeObjectSim(species = species)

        cl <- sim.dt@cl@cl[sim.dt@cl@cl$taxon == species,]
        cell.time <- 1:length(table(cl$time))
        cell.space <- 1:length(table(cl$space))
        cell.technical <- 1:length(table(cl$technical))
        cells <- expand.grid(cell.space, cell.technical, cell.time)

        cell.time.nms <-  names(table(cl$time))
        cell.space.nms <-  names(table(cl$space))
        cell.technical.nms <-  names(table(cl$technical))

        slots <- list(lenStruc = c('estim', 'rep'), lenVar = NA,
                      ageStruc = c('estim', 'rep'), ageVar = NA,
                      totalN   = c('estim', 'rep'), totalN = NA,
                      totalW   = c('estim', 'rep'), totalW = NA)



       for(i in 1:length(sim.dt@samples)){

       cat('--------------------------------------------------\n')
       cat('Sample: ', i,'\n')
       cat('--------------------------------------------------\n')

            # The bayesian fit for sample i
            fit <-run.bayes(sim.dt@samples[[i]], sim.dt@cl,fit = NULL, do.predict = FALSE, species,
                         season.definition = sim.dt@setup.args$season.definition,
                         burnin = burnin, thin = thin, nmcmc = nmcmc,
                         ageMin= sim.dt@setup.args$ageMin, ageMax = sim.dt@setup.args$ageMax,
                         wgl.eqn = sim.dt@setup.args$wgl.eqn,
                         age.covariates = sim.dt@setup.args$age.covariates,
                         weight.covariates = sim.dt@setup.args$weight.covariates,
                         prediction.cells = list(areas=NULL,gears=NULL,seasons= NULL),
                         arealist = NULL,
                         length.list = sim.dt@setup.args$length.list,
                         adj = NULL)

            # The prediction among cells.
            for(j in 1:dim(cells)[1]){

                pred.cell <- predict.by.clDataCons(sim.dt@samples[[i]],sim.dt@cl, fit=fit$fit, species,
                         season.definition = sim.dt@setup.args$season.definition,
                         ageMin= sim.dt@setup.args$ageMin, ageMax = sim.dt@setup.args$ageMax,
                         wgl.eqn = sim.dt@setup.args$wgl.eqn,
                         age.covariates = sim.dt@setup.args$age.covariates,
                         weight.covariates = sim.dt@setup.args$weight.covariates,
                         prediction.cells = list(areas = cells[j,1] ,gears=cells[j,2], seasons= cells[j,3]),
                         arealist = NULL,
                         length.list = sim.dt@setup.args$length.list)

               cat('       - - - - - - - - \n')
               cat( 'time: ', cells[j,3], ' space: ', cells[j,1], ' technical: ', cells[j,2], '\n')
               cat('       - - - - - - - - \n')
                res.lan <- pred.cell[['dbeObject.list']][['dbeTotcatch.land']]
                lan.wga <- pred.cell[['dbeObject.list']][['dbeMean.wga.land']]

                res.lan@lenStruc$estim$time <- cell.time.nms[cells[j,3]]
                res.lan@lenStruc$estim$space <- cell.space.nms[cells[j,1]]
                res.lan@lenStruc$estim$technical <- cell.technical.nms[cells[j,2]]
                res.lan@lenStruc$rep$time <- cell.time.nms[cells[j,3]]
                res.lan@lenStruc$rep$space <- cell.space.nms[cells[j,1]]
                res.lan@lenStruc$rep$technical <- cell.technical.nms[cells[j,2]]
                res.lan@lenVar$time <- cell.time.nms[cells[j,3]]
                res.lan@lenVar$space <- cell.space.nms[cells[j,1]]
                res.lan@lenVar$technical <- cell.technical.nms[cells[j,2]]
                res.lan@ageStruc$estim$time <- cell.time.nms[cells[j,3]]
                res.lan@ageStruc$estim$space <- cell.space.nms[cells[j,1]]
                res.lan@ageStruc$estim$technical <- cell.technical.nms[cells[j,2]]
                res.lan@ageStruc$rep$time <- cell.time.nms[cells[j,3]]
                res.lan@ageStruc$rep$space <- cell.space.nms[cells[j,1]]
                res.lan@ageStruc$rep$technical <- cell.technical.nms[cells[j,2]]
                res.lan@ageVar$time <- cell.time.nms[cells[j,3]]
                res.lan@ageVar$space <- cell.space.nms[cells[j,1]]
                res.lan@ageVar$technical <- cell.technical.nms[cells[j,2]]

                # Total N.
                totalNr <- aggregate(list( 'value' = res.lan@lenStruc$rep$value),
                                    list(time = res.lan@lenStruc$rep$time, space = res.lan@lenStruc$rep$space,
                                        technical = res.lan@lenStruc$rep$technical, iter = res.lan@lenStruc$rep$iter), sum)
                res.lan@totalN$rep <- totalNr[, c('time', 'space', 'technical', 'value', 'iter')]
                totalNe <- aggregate(list('value' = totalNr$value), list(time = totalNr$time, space = totalNr$space, technical = totalNr$technical), median)
                res.lan@totalN$estim <- totalNe
                totalNv <- aggregate(list('value' = totalNr$value), list(time = totalNr$time, space = totalNr$space, technical = totalNr$technical), var)
                res.lan@totalNvar <- totalNv

                # Total W.
                totalWr <- aggregate(list('value' = res.lan@ageStruc$rep$value*lan.wga@ageStruc$rep$value),
                                list(time = res.lan@ageStruc$rep$time, space = res.lan@ageStruc$rep$space,
                                technical = res.lan@ageStruc$rep$technical, iter = res.lan@ageStruc$rep$iter), sum)
                totalWr$value <- totalWr$value/1000
                res.lan@totalW$rep <- totalWr[, c('time', 'space', 'technical', 'value', 'iter')]
                totalWe <- aggregate(list('value' = totalWr$value), list(time = totalWr$time, space = totalWr$space, technical = totalWr$technical), median)
                res.lan@totalW$estim <- totalWe
                totalWv <- aggregate(list('value' = totalWr$value), list(time = totalWr$time, space = totalWr$space, technical = totalWr$technical), var)
                res.lan@totalWvar <- totalWv

                # aggregate to the total.


                for(s in 1:length(slots)){

                    x <- slot(res@LAN, names(slots)[s])


                    if(class(x) != 'list'){
                        ind <- 1:dim(slot(res@LAN, names(slots)[s]))[1]
                        if(i == 1) ind <- -1
                        slot(res@LAN, names(slots)[s]) <- rbind(slot(res@LAN, names(slots)[s])[ind,], cbind(sample = i, slot(res.lan, names(slots)[s])))
                    }
                    else{
                        ind <- 1:dim(slot(res@LAN, names(slots)[s])$estim)[1]
                        if(i == 1) ind <- -1
                        slot(res@LAN, names(slots)[s])$estim <- rbind(slot(res@LAN, names(slots)[s])$estim[ind,], cbind(sample = i, slot(res.lan, names(slots)[s])$estim))
                        if(rep == TRUE){
                            ind <- 1:dim(slot(res@LAN, names(slots)[s])$estim)[1]
                            if(i == 1) ind <- -1
                                slot(res@LAN, names(slots)[s])$rep <- rbind(slot(res@LAN, names(slots)[s])$rep[ind,], cbind(sample = i, slot(res.lan, names(slots)[s])$rep))
                       }
                   }
             }
         }
    }
    return(res)
}