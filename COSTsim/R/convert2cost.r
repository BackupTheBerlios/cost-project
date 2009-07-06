convert2cost<-function(sim.object,gear.list,area.list,species, year){
if(is.null(gear.list))gear.list<-1
if(is.null(area.list))area.list<-1
  months.num <- 1:12
names(months.num ) <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
sim.obs<-sim.ml<-NULL
  if(length(sim.object$obs@tr$trpCode)>0){
# OBS: TR ---------------------------------------------------------------------------
ntr <- length(sim.object$obs@tr$trpCode)
obs.tr <- data.frame(sampType   = rep('S', ntr),
                     landCtry   = rep('XYZ', ntr),
                     vslFlgCtry = rep('XYZ', ntr),
                     year       = rep(year, ntr),                
                     proj       = rep('NK', ntr),
                     trpCode    = sim.object$obs@tr$trpCode,
                     vslLen     = rep(NA, ntr),
                     vslPwr     = rep(NA, ntr),
                     vslSize    = rep(NA, ntr),
                     vslType    = rep(NA, ntr),
                     harbour    = rep(NA, ntr),
                     foNum      = sim.object$obs@tr$foNum,
                     daysAtSea  = rep(NA, ntr),
                     vslId      = rep(NA, ntr),
                     sampCtry   = rep(NA, ntr),
                     sampMeth   = 'observer')


# OBS: HH ---------------------------------------------------------------------------
nhh <- length(sim.object$obs@hh$trpCode)
obs.hh <- data.frame(sampType   = rep('S', nhh),
                     landCtry   = rep('XYZ', nhh),
                     vslFlgCtry = rep('XYZ', nhh),
                     year       = rep(year, nhh),
                     proj       = rep('NK', nhh),
                     trpCode    = sim.object$obs@hh$trpCode,
                     staNum     = sim.object$obs@hh$staNum,
                     foVal      = rep('V', nhh),
                     aggLev     = rep('H', nhh),
                     catReg     = 'All',
                     sppReg     = 'Par',
                     date       = sim.object$obs@hh$date,
                     time       = rep(NA, nhh),
                     foDur      = rep(NA, nhh),
                     latIni     = rep(NA, nhh),
                     lonIni     = rep(NA, nhh),
                     latFin     = rep(NA, nhh),
                     lonFin     = rep(NA, nhh),
                     area       = area.list[sim.object$obs@hh$area],
                     rect       = rep(NA, nhh),
                     subRect    = rep(NA, nhh),
                     foDep      = rep(NA, nhh),
                     waterDep   = rep(NA, nhh),
                     foCatNat   = rep(NA, nhh),
                     foCatEu5   = gear.list[sim.object$obs@hh$foCatEu5],
                     foCatEu6   = rep(NA, nhh),
                     meshSize   = rep(NA, nhh),
                     selDev     = rep('0', nhh),
                     meshSizeSelDev =  rep(NA, nhh))

# OBS: SL ---------------------------------------------------------------------------
#a <- sim.object$obs@sl

#b <- aggregate(list(wt = a$wt, subSampWt = a$subSampWt),
#    list(catchCat = a$catchCat, spp = a$spp, trpCode = a$trpCode,
#    staNum = ifelse(is.na(a$staNum), 'kk', a$staNum),
#    commCat = ifelse(is.na(a$commCat), 'kk', a$commCat)), sum)


nsl <- length(sim.object$obs@sl$trpCode)
obs.sl <- data.frame(sampType   = rep('S', nsl),
                     landCtry   = rep('XYZ', nsl),
                     vslFlgCtry = rep('XYZ', nsl),
                     year       = rep(year, nsl),
                     proj       = rep('NK', nsl),
                     trpCode    = sim.object$obs@sl$trpCode,
                     staNum     = sim.object$obs@sl$staNum,
                     spp        = rep(species, nsl),
                     catchCat   = sim.object$obs@sl$catchCat,
                     landCat    = rep(NA, nsl),
                     commCatScl = rep(NA, nsl),
                     commCat    = as.character(sim.object$obs@sl$commCat),
                     subSampCa  = rep(NA, nsl),
                     sex        = rep(NA, nsl),
                     wt         = sim.object$obs@sl$wt,
                     subSampWt  = sim.object$obs@sl$subSampWt,
                     lenCode    = rep(NA, nsl))



# OBS: HL ---------------------------------------------------------------------------
nhl <- length(sim.object$obs@hl$trpCode)
obs.hl <- data.frame(sampType   = rep('S', nhl),
                     landCtry   = rep('XYZ', nhl),
                     vslFlgCtry = rep('XYZ', nhl),
                     year       = rep(year, nhl),
                     proj       = rep('NK', nhl),
                     trpCode    = sim.object$obs@hl$trpCode,
                     staNum     = sim.object$obs@hl$staNum,
                     spp        = rep(species, nhl),
                     catchCat   = sim.object$obs@hl$catchCat,
                     landCat    = rep(NA, nhl),
                     commCatScl = rep(NA, nhl),
                     commCat    = as.character(sim.object$obs@hl$commCat),
                     subSampCa  = rep(NA, nhl),
                     sex        = rep(NA, nhl),
                     lenCls     = sim.object$obs@hl$lenCls,
                     lenNum     = sim.object$obs@hl$lenNum)

# OBS: CA ---------------------------------------------------------------------------
nca <- length(sim.object$obs@ca$trpCode)

# date
ca.date <- as.Date(rep(NA, nca))
for(i in 1:nca){
    ca.date[i] <- as.Date(sim.object$obs@hh[sim.object$obs@hh$trpCode == sim.object$obs@ca$trpCode[i],]$date[1])
}

obs.ca <- data.frame(sampType   = rep('S', nca),
                     landCtry   = rep('XYZ', nca),
                     vslFlgCtry = rep('XYZ', nca),
                     year       = rep(year, nca),
                     proj       = rep('NK', nca),
                     trpCode    = sim.object$obs@ca$trpCode,
                     staNum     = rep(NA, nca),
                     quarter    = as.numeric(substr(quarters(ca.date),2,2)),
                     month      = months.num[months(ca.date, abb = TRUE)],
                     spp        = rep(species, nca),
                     sex        = rep(NA, nca),
                     catchCat   = sim.object$obs@ca$catchCat,
                     landCat    = rep(NA, nca),
                     commCatScl = rep(NA, nca),
                     commCat    = rep(NA, nca),
                     stock      = rep('her-vian', nca),
                     area       = rep(NA, nca),
                     rect       = rep(NA, nca),
                     subRect    = rep(NA, nca),
                     lenCls     = sim.object$obs@ca$lenCls,
                     age        = sim.object$obs@ca$age,
                     fishId     = 1:nca,
                     lenCode    = rep(NA, nca),
                     ageMeth    = rep(NA, nca),
                     plusGrp    = rep(NA, nca),
                     otoWt      = rep(NA, nca),
                     otoSide    = rep(NA, nca),
                     indWt      = sim.object$obs@ca$wt,
                     matMeth    = rep(NA, nca),
                     matScale   = rep(NA, nca),
                     matStage   = rep(NA, nca))


# OBS: HH ---------------------------------------------------------------------------
# sl observations which are taken at trip and not metier level.
hh0 <- unique(obs.sl[which(is.na(obs.sl$staNum)),1:7])
nhh <- dim(hh0)[1]

date0 <- numeric(nhh)
area0 <- numeric(nhh)
foC0 <- numeric(nhh)
for(i in 1:nhh){
    date0[i] <- as.character(obs.hh[obs.hh$trpCode == hh0[i, 'trpCode'],][1,'date'])
    area0[i] <- as.character(obs.hh[obs.hh$trpCode == hh0[i, 'trpCode'],][1,'area'])
    foC0[i] <- as.character(obs.hh[obs.hh$trpCode == hh0[i, 'trpCode'],][1,'foCatEu5'])
}


obs.hh0 <- data.frame(sampType   = rep('S', nhh),
                     landCtry   = rep('XYZ', nhh),
                     vslFlgCtry = rep('XYZ', nhh),
                     year       = rep(year, nhh),
                     proj       = rep('NK', nhh),
                     trpCode    = hh0$trpCode,
                     staNum     = hh0$staNum,
                     foVal      = rep('V', nhh),
                     aggLev     = rep('T', nhh),
                     catReg     = 'All',
                     sppReg     = 'Par',
                     date       = date0,
                     time       = rep(NA, nhh),
                     foDur      = rep(NA, nhh),
                     latIni     = rep(NA, nhh),
                     lonIni     = rep(NA, nhh),
                     latFin     = rep(NA, nhh),
                     lonFin     = rep(NA, nhh),
                     area       = area0,
                     rect       = rep(NA, nhh),
                     subRect    = rep(NA, nhh),
                     foDep      = rep(NA, nhh),
                     waterDep   = rep(NA, nhh),
                     foCatNat   = rep(NA, nhh),
                     foCatEu5   = foC0,
                     foCatEu6   = rep(NA, nhh),
                     meshSize   = rep(NA, nhh),
                     selDev     = rep('0', nhh),
                     meshSizeSelDev =  rep(NA, nhh))

# obs object
sim.obs <- csData(tr = obs.tr, hh = rbind(obs.hh, obs.hh0), hl = obs.hl, sl = obs.sl, ca = obs.ca,desc = 'Simulated data')
validObject(sim.obs)
}


  
#-------------------------------------------------------------------------------
# ML
#-------------------------------------------------------------------------------
  if(length(sim.object$ml@tr$trpCode)>0){

# ML: TR ---------------------------------------------------------------------------
ntr <- length(sim.object$ml@tr$trpCode)
ml.tr <- data.frame(sampType   = rep('M', ntr),
                     landCtry   = rep('XYZ', ntr),
                     vslFlgCtry = rep('XYZ', ntr),
                     year       = rep(year, ntr),
                     proj       = rep('NK', ntr),
                     trpCode    = sim.object$ml@tr$trpCode,
                     vslLen     = rep(NA, ntr),
                     vslPwr     = rep(NA, ntr),
                     vslSize    = rep(NA, ntr),
                     vslType    = rep(NA, ntr),
                     harbour    = rep(NA, ntr),
                     foNum      = sim.object$ml@tr$foNum,
                     daysAtSea  = rep(NA, ntr),
                     vslId      = rep(NA, ntr),
                     sampCtry   = rep(NA, ntr),
                     sampMeth   = 'observer')


# ML: HH ---------------------------------------------------------------------------
nhh <- length(sim.object$ml@hh$trpCode)
ml.hh <- data.frame(sampType   = rep('M', nhh),
                     landCtry   = rep('XYZ', nhh),
                     vslFlgCtry = rep('XYZ', nhh),
                     year       = rep(year, nhh),
                     proj       = rep('NK', nhh),
                     trpCode    = sim.object$ml@hh$trpCode,
                     staNum     = sim.object$ml@hh$staNum,
                     foVal      = rep('V', nhh),
                     aggLev     = rep('H', nhh),
                     catReg     = 'All',
                     sppReg     = 'Par',
                     date       = sim.object$ml@hh$date,
                     time       = rep(NA, nhh),
                     foDur      = rep(NA, nhh),
                     latIni     = rep(NA, nhh),
                     lonIni     = rep(NA, nhh),
                     latFin     = rep(NA, nhh),
                     lonFin     = rep(NA, nhh),
                     area       = area.list[sim.object$ml@hh$area],
                     rect       = rep(NA, nhh),
                     subRect    = rep(NA, nhh),
                     foDep      = rep(NA, nhh),
                     waterDep   = rep(NA, nhh),
                     foCatNat   = rep(NA, nhh),
                     foCatEu5   = gear.list[sim.object$ml@hh$foCatEu5],
                     foCatEu6   = rep(NA, nhh),
                     meshSize   = rep(NA, nhh),
                     selDev     = rep('0', nhh),
                     meshSizeSelDev =  rep(NA, nhh))

# ML: SL ---------------------------------------------------------------------------
a <- sim.object$ml@sl

b <- aggregate(list(wt = a$wt, subSampWt = a$subSampWt),
    list(catchCat = a$catchCat, spp = a$spp, trpCode = a$trpCode,
    commCat = as.character(a$commCat)), sum)


nsl <- dim(b)[1]
ml.sl <- data.frame(sampType   = rep('M', nsl),
                     landCtry   = rep('XYZ', nsl),
                     vslFlgCtry = rep('XYZ', nsl),
                     year       = rep(year, nsl),
                     proj       = rep('NK', nsl),
                     trpCode    = b$trpCode,
                     staNum     = rep(NA, nsl), #sim.object$ml@sl$staNum,
                     spp        = rep(species, nsl),
                     catchCat   = b$catchCat,
                     landCat    = rep(NA, nsl),
                     commCatScl = rep(NA, nsl),
                     commCat    = b$commCat,
                     subSampCa  = rep(NA, nsl),
                     sex        = rep(NA, nsl),
                     wt         = b$wt,
                     subSampWt  = b$subSampWt,
                     lenCode    = rep(NA, nsl))



# ML: HL ---------------------------------------------------------------------------
a <- sim.object$ml@hl
b <- aggregate(list(lenNum = a$lenNum),
    list(catchCat = a$catchCat, trpCode = a$trpCode, commCat = as.character(a$commCat), lenCls = a$lenCls), sum)

nhl <- dim(b)[1]
ml.hl <- data.frame(sampType   = rep('M', nhl),
                     landCtry   = rep('XYZ', nhl),
                     vslFlgCtry = rep('XYZ', nhl),
                     year       = rep(year, nhl),
                     proj       = rep('NK', nhl),
                     trpCode    = b$trpCode,
                     staNum     = rep(NA, nhl),
                     spp        = rep(species, nhl),
                     catchCat   = b$catchCat,
                     landCat    = rep(NA, nhl),
                     commCatScl = rep(NA, nhl),
                     commCat    = b$commCat,
                     subSampCa  = rep(NA, nhl),
                     sex        = rep(NA, nhl),
                     lenCls     = b$lenCls,
                     lenNum     = b$lenNum)

# ML: CA ---------------------------------------------------------------------------
nca <- length(sim.object$ml@ca$trpCode)

# date
ca.date <- as.Date(rep(NA, nca))
for(i in 1:nca){
    ca.date[i] <- as.Date(sim.object$ml@hh[sim.object$ml@hh$trpCode == sim.object$ml@ca$trpCode[i],]$date)
}

ml.ca <- data.frame(sampType   = rep('M', nca),
                     landCtry   = rep('XYZ', nca),
                     vslFlgCtry = rep('XYZ', nca),
                     year       = rep(year, nca),
                     proj       = rep('NK', nca),
                     trpCode    = sim.object$ml@ca$trpCode,
                     staNum     = rep(NA, nca),
                     quarter    = as.numeric(substr(quarters(ca.date),2,2)),
                     month      = months.num[months(ca.date, abb = TRUE)],
                     spp        = rep(species, nca),
                     sex        = rep(NA, nca),
                     catchCat   = sim.object$ml@ca$catchCat,
                     landCat    = rep(NA, nca),
                     commCatScl = rep(NA, nca),
                     commCat    = rep(NA, nca),
                     stock      = rep('her-vian', nca),
                     area       = rep(NA, nca),
                     rect       = rep(NA, nca),
                     subRect    = rep(NA, nca),
                     lenCls     = sim.object$ml@ca$lenCls,
                     age        = sim.object$ml@ca$age,
                     fishId     = 1:nca,
                     lenCode    = rep(NA, nca),
                     ageMeth    = rep(NA, nca),
                     plusGrp    = rep(NA, nca),
                     otoWt      = rep(NA, nca),
                     otoSide    = rep(NA, nca),
                     indWt      = sim.object$ml@ca$wt,
                     matMeth    = rep(NA, nca),
                     matScale   = rep(NA, nca),
                     matStage   = rep(NA, nca))


# ml object
sim.ml <- csData(tr = ml.tr, hh = ml.hh, hl = ml.hl, sl = ml.sl, ca = ml.ca, desc = 'Simulated data')
validObject(sim.ml)
}

sim.cs<-NULL
if(!is.null(sim.obs)&!is.null(sim.ml))sim.cs <- rbind2(sim.obs, sim.ml)
  else {if(!is.null(sim.obs))sim.cs<-sim.obs
        if(!is.null(sim.ml))sim.cs<-sim.ml}
        
sim.cs

}
