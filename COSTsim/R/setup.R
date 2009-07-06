setup<-function(data.args){
    new.setup<-function(seaslist,arealist,gearlist,species,cov.list,ntrip=50,ndisc=30,div,ageMin,ageMax)
    
    datalist <- data.args$datalist
    species  <- data.args$species
    cov.list <- data.args$cov.list
    ntrip    <- data.args$ntrip
    ndisc    <- data.args$ndisc
    div      <- data.args$div
    ageMin   <- data.args$ageMin
    ageMax   <- data.args$ageMax

    res <- new.setup(seaslist,arealist,gearlist,species,cov.list,ntrip=50,ndisc=30,div,ageMin,ageMax)
    return(res)
}
