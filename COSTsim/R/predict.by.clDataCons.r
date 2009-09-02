#-------------------------------------------------------------------------------
# Dorleta Garcia   - 27/08/2009 08:25:48
# 
#
# predict.by.clDataCons: Based on the predict part of the run.bayes function and
#   adapted to work with 'consolidated' data sets.
# mbe2dbe.dor: modification of the mbe2dbe function there was a problem in the original.
#-------------------------------------------------------------------------------


predict.by.clDataCons <- function (csdata, cldata, fit, species,
    season.definition = "quarter",
    ageMin = 0, ageMax = 20, wgl.eqn = NULL, age.covariates = NULL,
    weight.covariates = NULL, prediction.cells = list(areas = "ALL",
        gears = "ALL", seasons = "ALL"), arealist = NULL, length.list)
        
{
    obsdata <- mldata <- NULL
    if (sum(csdata@tr$sampType == "M") > 0)
        mldata <- COSTcore:::subset(csdata, subset = (sampType == "M"))
    if (sum(csdata@tr$sampType == "S") > 0)
        obsdata <- COSTcore:::subset(csdata, subset = (sampType == "S"))
    if (is.null(mldata) & (is.null(obsdata))) {
        print("no data of type M or S")
        return
    }
    nseas <- 0
    if (season.definition == "quarter")
        nseas <- 4
    if (season.definition == "month")
        nseas <- 12
    if (nseas == 0 & age.covariates$season) {
        print("invalid time strata")
        return
    }
    models <- make.reg.models(age.covariates, weight.covariates)
    l.int <- log(seq(length.list$minl, length.list$maxl, length.list$int))
    if (!is.null(wgl.eqn)) {
        if (!is.null(wgl.eqn$int$seas)) {
            if (nseas == 4 & length(wgl.eqn$int$seas) == 12)
                wgl.eqn$int$seas <- c(mean(wgl.eqn$int$seas[1:3]),
                  mean(wgl.eqn$int$seas[4:6]), mean(wgl.eqn$int$seas[7:9]),
                  mean(wgl.eqn$int$seas[10:12]))
        }
    }
    input.data <- read.cost.data(obsdata, mldata, species, wgl.eqn,
        season.definition)
    trip.alk <- rep(1:input.data$n_trip_mland, input.data$num_alk_mland)
    if (!age.covariates$season) {
        input.data$season_mland <- rep(1, length(input.data$season_mland))
        input.data$seas <- rep(1, length(input.data$seas))
    }
    trip.l1 <- rep(1:input.data$n_trip_mland, input.data$num_trip_mland)
    trip.l <- rep(trip.l1, input.data$num_size_mland)
    input.data$trip.l <- trip.l
    input.data$trip.alk <- trip.alk
    if (!is.null(input.data$alk_a_disc)) {
        input.data$alk_a_disc[input.data$alk_a_disc < ageMin] <- ageMin
        input.data$alk_a_disc[input.data$alk_a_disc > ageMax] <- ageMax
    }
    if (!is.null(input.data$alk_a_mland)) {
        input.data$alk_a_mland[input.data$alk_a_mland < ageMin] <- ageMin
        input.data$alk_a_mland[input.data$alk_a_mland > ageMax] <- ageMax
    }
    if (!is.null(input.data$sampsize_disc))
        input.data$sampsize_disc[is.na(input.data$sampsize_disc)] <- 0
    if (!is.null(input.data$haulsize_disc))
        input.data$haulsize_disc[is.na(input.data$haulsize_disc)] <- 0
    input.data$cens.mu = c(1, 30, 3.4)
    input.data$cens.tau = c(1, 1, 1)
    input.data$cens.pri = c(3.4, 0.001, 1e-04, 0.001)
    seas.code <- sim.codecov(input.data$seas)
    gear.code <- sim.codecov(input.data$gear)
    area.code <- sim.codecov(input.data$area, arealist)



    predict <- dbeObject.list <- NULL

        seas.pred <- sim.codecov(prediction.cells$seas, seas.code$covlist)
        gear.pred <- sim.codecov(prediction.cells$gear, gear.code$covlist)
        area.pred <- sim.codecov(prediction.cells$area, area.code$covlist)
        if (length(prediction.cells$seas) == 1)
            if (prediction.cells$seas == "ALL")
                prediction.cells$seas <- seas.code$covlist
        if (length(prediction.cells$gear) == 1)
            if (prediction.cells$gear == "ALL")
                prediction.cells$gear <- gear.code$covlist
        if (length(prediction.cells$area) == 1)
            if (prediction.cells$area == "ALL")
                prediction.cells$area <- area.code$covlist
        newpred <- expand.grid(1, prediction.cells$seas, prediction.cells$gear,
            prediction.cells$area)
        tseas <- code.cov(newpred[, 2], seas.code$covlist)
        tgear <- code.cov(newpred[, 3], gear.code$covlist)
        tarea <- code.cov(newpred[, 4], area.code$covlist)
        tyear <- rep(1, length(tseas))
        if (!models$agemodel$Int$area)
            tarea <- rep(1, length(tseas))
        if (!models$agemodel$Int$seas)
            tseas <- rep(1, length(tseas))
        if (!models$agemodel$Int$gear)
            tgear <- rep(1, length(tseas))
        landings <- NULL
        land.seas <- cldata@cl$time
        
 land.seas <- as.numeric(matrix(unlist(strsplit(as.character(land.seas), ' - ')),2)[2,])
 
  #      if (season.definition == "quarter")
  #          land.seas <- cldata@cl$quarter
  #      if (season.definition == "month")
  #          land.seas <- cldata@cl$month
        if (season.definition != "quarter" & season.definition !=
            "month" & age.covariates$season) {
            print("invalid season.definition")
            return
        }
        for (i in 1:length(newpred[, 1])) {
            ind <- cldata@cl$taxon == species
            if (models$agemodel$Int$seas)
                ind <- ind & (land.seas == newpred[i, 2])
            if (models$agemodel$Int$gear)
                ind <- ind & (cldata@cl$technical == newpred[i,
                  3])
            if (models$agemodel$Int$area)
                ind <- ind & (cldata@cl$space == newpred[i, 4])
            landings <- c(landings, sum(cldata@cl$landWt[ind]))
        }
        landings <- 1000 * landings
        predict <- predict.fit.COST(fit, fit$COST.list, tyear,
            tseas, tgear, tarea, landings, tseas, t2.year = NULL,
            t2.seas = NULL, t2.gear = NULL, t2.area = NULL, burnin = 0,
            nMC = 100, l.int = l.int, par.haulsize = NULL)
        dbeObject.list <- mbe2dbe.dor(predict, species)

    list(predict = predict, dbeObject.list = dbeObject.list)
}


mbe2dbe.dor <- function (mbeoutput, species = as.character(NA)) 
{
    if (nrow(mbeoutput$cov) == 1) {
        time <- space <- technical <- "all"
    }
    else {
        time <- paste(mbeoutput$cov$year, mbeoutput$cov$seas, 
            sep = " - ")
        space <- mbeoutput$cov$area
        technical <- mbeoutput$cov$gear
    }
    dimnames(mbeoutput$totcatch.land) <- list(round(exp(mbeoutput$l.int)), 
        mbeoutput$avec, 1:(dim(mbeoutput$totcatch.land)[3]))
    dimnames(mbeoutput$mean.lga.land) <- dimnames(mbeoutput$mean.wga.land) <- list(mbeoutput$avec, 
        1:(dim(mbeoutput$totcatch.land)[3]))
    if (!is.null(mbeoutput$totcatch.disc)) {
        dimnames(mbeoutput$totcatch.disc) <- dimnames(mbeoutput$totcatch.land)
        dimnames(mbeoutput$mean.lga.disc) <- dimnames(mbeoutput$mean.wga.disc) <- list(mbeoutput$avec, 
            1:(dim(mbeoutput$totcatch.land)[3]))
    }
    dbeTotcatch.land <- dbeTotcatch.disc <- dbeMean.lga.land <- dbeMean.wga.land <- dbeMean.lga.disc <- dbeMean.wga.disc <- new("dbeOutput", 
        species = species)
    dbeTotcatch.land@param <- dbeTotcatch.disc@param <- "N"
    dbeMean.lga.land@param <- dbeMean.lga.disc@param <- "length"
    dbeMean.wga.land@param <- dbeMean.wga.disc@param <- "weight"
    dbeTotcatch.land@catchCat <- dbeMean.lga.land@catchCat <- dbeMean.wga.land@catchCat <- "LAN"
    dbeTotcatch.disc@catchCat <- dbeMean.lga.disc@catchCat <- dbeMean.wga.disc@catchCat <- "DIS"
    atLgth <- apply(mbeoutput$totcatch.land, c(1, 3), sum, na.rm = TRUE)
    atAge <- apply(mbeoutput$totcatch.land, 2:3, sum, na.rm = TRUE)
    expLgth <- expand.grid(dimnames(atLgth))
    expAge <- expand.grid(dimnames(atAge))
    dbeTotcatch.land@lenStruc$rep <- data.frame(cbind(time = time, 
        space = space, technical = technical, length = as.character(expLgth[, 
            1])), value = as.vector(atLgth), iter = expLgth[, 
        2])
    dbeTotcatch.land@ageStruc$rep <- data.frame(cbind(time = time, 
        space = space, technical = technical, age = as.character(expAge[, 
            1])), value = as.vector(atAge), iter = expAge[, 2])
    dbeTotcatch.land = fillTabFromRep(dbeTotcatch.land)
    dbeTotcatch.disc <- NULL
    if (!is.null(mbeoutput$totcatch.disc)) {
        atLgth <- apply(mbeoutput$totcatch.disc, c(1, 3), sum, 
            na.rm = TRUE)
        atAge <- apply(mbeoutput$totcatch.disc, 2:3, sum, na.rm = TRUE)
        expLgth <- expand.grid(dimnames(atLgth))
        expAge <- expand.grid(dimnames(atAge))
        dbeTotcatch.disc@lenStruc$rep <- data.frame(cbind(time = time, 
            space = space, technical = technical, length = as.character(expLgth[, 
                1])), value = as.vector(atLgth), iter = expLgth[, 
            2])
        dbeTotcatch.disc@ageStruc$rep <- data.frame(cbind(time = time, 
            space = space, technical = technical, age = as.character(expAge[, 
                1])), value = as.vector(atAge), iter = expAge[, 
            2])
        dbeTotcatch.disc <- fillTabFromRep(dbeTotcatch.disc)
    }
    expAge <- expand.grid(dimnames(mbeoutput$mean.lga.land))
    dbeMean.lga.land@ageStruc$rep <- data.frame(cbind(time = time, 
        space = space, technical = technical, age = as.character(expAge[, 
            1])), value = as.vector(mbeoutput$mean.lga.land), 
        iter = expAge[, 2])
    dbeMean.lga.land = fillTabFromRep(dbeMean.lga.land)
    expAge <- expand.grid(dimnames(mbeoutput$mean.wga.land))
    dbeMean.wga.land@ageStruc$rep <- data.frame(cbind(time = time, 
        space = space, technical = technical, age = as.character(expAge[, 
            1])), value = as.vector(mbeoutput$mean.wga.land), 
        iter = expAge[, 2])
    dbeMean.wga.land = fillTabFromRep(dbeMean.wga.land)
    dbeMean.lga.disc <- NULL
    if (!is.null(mbeoutput$mean.lga.disc)) {
        expAge <- expand.grid(dimnames(mbeoutput$mean.lga.disc))
        dbeMean.lga.disc@ageStruc$rep <- data.frame(cbind(time = time, 
            space = space, technical = technical, age = as.character(expAge[, 
                1])), value = as.vector(mbeoutput$mean.lga.disc), 
            iter = expAge[, 2])
        dbeMean.lga.disc = fillTabFromRep(dbeMean.lga.disc)
    }
    dbeMean.wga.disc <- NULL
    if (!is.null(mbeoutput$mean.wga.disc)) {
        expAge <- expand.grid(dimnames(mbeoutput$mean.wga.disc))
        dbeMean.wga.disc@ageStruc$rep <- data.frame(cbind(time = time, 
            space = space, technical = technical, age = as.character(expAge[, 
                1])), value = as.vector(mbeoutput$mean.wga.disc), 
            iter = expAge[, 2])
        dbeMean.wga.disc = fillTabFromRep(dbeMean.wga.disc)
    }
    return(list(dbeTotcatch.land = dbeTotcatch.land, dbeTotcatch.disc = dbeTotcatch.disc, 
        dbeMean.lga.land = dbeMean.lga.land, dbeMean.wga.land = dbeMean.wga.land, 
        dbeMean.lga.disc = dbeMean.lga.disc, dbeMean.wga.disc = dbeMean.wga.disc))
}
