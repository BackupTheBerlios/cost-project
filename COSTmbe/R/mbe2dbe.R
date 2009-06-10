mbe2dbe <- function(mbeoutput,species=as.character(NA)) {

if (nrow(mbeoutput$cov)>1) {
time <- space <- technical <- "all"
} else {
#stratification definition
time <- paste(mbeoutput$cov$year,mbeoutput$cov$seas,sep=" - ")
space <- mbeoutput$cov$area
technical <- mbeoutput$cov$gear
}

#mbeoutput's dimnames
dimnames(mbeoutput$totcatch.land) <- dimnames(mbeoutput$totcatch.disc) <- list(round(exp(mbeoutput$l.int)),mbeoutput$avec,1:(dim(mbeoutput$totcatch.land)[3]))
dimnames(mbeoutput$mean.lga.land) <- dimnames(mbeoutput$mean.wga.land) <- list(mbeoutput$avec,1:(dim(mbeoutput$totcatch.land)[3]))
dimnames(mbeoutput$mean.lga.disc) <- dimnames(mbeoutput$mean.wga.disc) <- list(mbeoutput$avec,1:(dim(mbeoutput$totcatch.land)[3]))

#dbeOutput objects are created
dbeTotcatch.land <- dbeTotcatch.disc <- dbeMean.lga.land <- dbeMean.wga.land <- dbeMean.lga.disc <- dbeMean.wga.disc <-new("dbeOutput",species=species)

#'param' slot is filled
dbeTotcatch.land@param <- dbeTotcatch.disc@param <- "N"
dbeMean.lga.land@param <- dbeMean.lga.disc@param <- "length"
dbeMean.wga.land@param <- dbeMean.wga.disc@param <- "weight"

#'catchCat' slot is filled
dbeTotcatch.land@catchCat <- dbeMean.lga.land@catchCat <- dbeMean.wga.land@catchCat <- "LAN"
dbeTotcatch.disc@catchCat <- dbeMean.lga.disc@catchCat <- dbeMean.wga.disc@catchCat <- "DIS"

#totcatch.land
atLgth <- apply(mbeoutput$totcatch.land,c(1,3),sum,na.rm=TRUE)
atAge <- apply(mbeoutput$totcatch.land,2:3,sum,na.rm=TRUE)
expLgth <- expand.grid(dimnames(atLgth)) ; expAge <- expand.grid(dimnames(atAge))

dbeTotcatch.land@lenStruc$rep <- data.frame(cbind(time=time,space=space,technical=technical,length=as.character(expLgth[,1]),value=as.vector(atLgth),iter=expLgth[,2]))
dbeTotcatch.land@ageStruc$rep <- data.frame(cbind(time=time,space=space,technical=technical,age=as.character(expAge[,1]),value=as.vector(atAge),iter=expAge[,2]))

#totcatch.disc
atLgth <- apply(mbeoutput$totcatch.disc,c(1,3),sum,na.rm=TRUE)
atAge <- apply(mbeoutput$totcatch.disc,2:3,sum,na.rm=TRUE)
expLgth <- expand.grid(dimnames(atLgth)) ; expAge <- expand.grid(dimnames(atAge))
dbeTotcatch.disc@lenStruc$rep <- data.frame(cbind(time=time,space=space,technical=technical,length=as.character(expLgth[,1]),value=as.vector(atLgth),iter=expLgth[,2]))
dbeTotcatch.disc@ageStruc$rep <- data.frame(cbind(time=time,space=space,technical=technical,age=as.character(expAge[,1]),value=as.vector(atAge),iter=expAge[,2]))

#mean.lga.land
expAge <- expand.grid(dimnames(mbeoutput$mean.lga.land))
dbeMean.lga.land@ageStruc$rep <- data.frame(cbind(time=time,space=space,technical=technical,age=as.character(expAge[,1]),value=as.vector(mbeoutput$mean.lga.land),iter=expAge[,2]))

#mean.wga.land
expAge <- expand.grid(dimnames(mbeoutput$mean.wga.land))
dbeMean.wga.land@ageStruc$rep <- data.frame(cbind(time=time,space=space,technical=technical,age=as.character(expAge[,1]),value=as.vector(mbeoutput$mean.wga.land),iter=expAge[,2]))

#mean.lga.disc
expAge <- expand.grid(dimnames(mbeoutput$mean.lga.disc))
dbeMean.lga.disc@ageStruc$rep <- data.frame(cbind(time=time,space=space,technical=technical,age=as.character(expAge[,1]),value=as.vector(mbeoutput$mean.lga.disc),iter=expAge[,2]))

#mean.wga.disc
expAge <- expand.grid(dimnames(mbeoutput$mean.wga.disc))
dbeMean.wga.disc@ageStruc$rep <- data.frame(cbind(time=time,space=space,technical=technical,age=as.character(expAge[,1]),value=as.vector(mbeoutput$mean.wga.disc),iter=expAge[,2]))

return(list(dbeTotcatch.land=dbeTotcatch.land,dbeTotcatch.disc=dbeTotcatch.disc,dbeMean.lga.land=dbeMean.lga.land,
            dbeMean.wga.land=dbeMean.wga.land,dbeMean.lga.disc=dbeMean.lga.disc,dbeMean.wga.disc=dbeMean.wga.disc))
}

