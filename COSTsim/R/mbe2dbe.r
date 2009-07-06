mbe2dbe.dor <- function (mbeoutput, species = as.character(NA)) 
{
    if (nrow(mbeoutput$cov) > 1) {
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
  #  dbeTotcatch.disc <- NULL
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
#    dbeMean.lga.disc <- NULL
    if (!is.null(mbeoutput$mean.lga.disc)) {
        expAge <- expand.grid(dimnames(mbeoutput$mean.lga.disc))
        dbeMean.lga.disc@ageStruc$rep <- data.frame(cbind(time = time, 
            space = space, technical = technical, age = as.character(expAge[, 
                1])), value = as.vector(mbeoutput$mean.lga.disc), 
            iter = expAge[, 2])
        dbeMean.lga.disc = fillTabFromRep(dbeMean.lga.disc)
    }
#    dbeMean.wga.disc <- NULL
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
