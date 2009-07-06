
truePopulation <- function(COST.data, by.fit, setup.args, desc = 'True population'){

     seas  <- tolower(setup.args$seasons)
     area  <- tolower(setup.args$areas)
     gear  <- tolower(setup.args$gears)
     
     if(is.null(seas)) seas <- 'all'
     if(is.null(area)) area <- 'all'
     if(is.null(gear)) gear <- 'all'
     
     age   <- setup.args$ageMin:setup.args$ageMax
     
  #   predict.cov <- make.predict.cov(seas,gear,area)

     l.int  <- seq(setup.args$length.list$minl, setup.args$length.list$maxl, setup.args$length.list$int)
     
     laa <- daa <- cbind(expand.grid(age = age, technical = gear, space = area, time = paste(setup.args$year, seas, sep = " - "))[,4:1], value = NA)
     lal <- dal <- cbind(expand.grid(length = l.int, technical = gear, space = area, time = paste(setup.args$year, seas, sep = " - "))[,4:1], value = NA)
     ltw <- ltn <- dtw <- dtn <- cbind(expand.grid(technical = gear, space = area, time = paste(setup.args$year, seas, sep = " - "))[,3:1], value = NA)

   #  lands       <- ltw
   #  lands$value <- setup.args$landings

     
     
     cfn <- 1 #1000*5      # correction factor numbers
     cfw <- 1 # 0.001/5    # correction factor weight
                                
     for(s in seas){
        for(a in area){
            for(g in gear){

            print(paste('Season: ', s, ' Area: ', a, ' Gear: ', g, sep = ""))
            
               prediction.cells <- list(areas= s,gears=g,seasons=s)
            
               truePop <- predict.true.bayes(COST.data@cs, COST.data@cl, fit=by.fit$fit, setup.args$species, season.definition= setup.args$season.definition, 
                           ageMin = setup.args$ageMin, ageMax = setup.args$ageMax, 
                            age.covariates = setup.args$age.covariates, weight.covariates = setup.args$weight.covariates, 
                            prediction.cells = prediction.cells, arealist = setup.args$arealist, length.list = setup.args$length.list)

                true   <- getparams(by.fit$fit,truePop$predict,simno = 1)
                true$caa$wga.land <- ifelse(true$caa$wga.land == 'Inf', 0,  true$caa$wga.land)
                true$caa$wga.disc <- ifelse(true$caa$wga.disc == 'Inf', 0,  true$caa$wga.disc)
                
                s1 <- paste(setup.args$year, ifelse(is.na(truePop$seas), s, truePop$seas), sep = " - ")
                a1 <- ifelse(is.na(truePop$area), a, truePop$area)
                g1 <- ifelse(is.na(truePop$gear), g, truePop$gear)
                
                laa[laa$time == s1 & laa$space == a1 & laa$technical == g1,'value'] <- true$caa$CAA.land*cfn
                lal[lal$time == s1 & lal$space == a1 & lal$technical == g1,'value'] <- apply(truePop$predict$totcatch.land,c(1,3),sum)*cfn
                ltw[ltw$time == s1 & ltw$space == a1 & ltw$technical == g1,'value'] <- sum(true$caa$CAA.land*true$caa$wga.land*cfw*cfn, na.rm = TRUE )
                ltn[ltn$time == s1 & ltn$space == a1 & ltn$technical == g1,'value'] <- sum(true$caa$CAA.land*cfn, na.rm = TRUE)

                if(!is.null(true$caa$CAA.disc)){
                    daa[daa$time == s1 & daa$space == a1 & daa$technical == g1,'value'] <- true$caa$CAA.disc*cfn
                    dal[dal$time == s1 & dal$space == a1 & dal$technical == g1,'value'] <- apply(truePop$predict$totcatch.disc,c(1,3),sum)*cfn
                    dtw[dtw$time == s1 & dtw$space == a1 & dtw$technical == g1,'value'] <- sum(true$caa$CAA.disc*true$caa$wga.discw*cfn, na.rm = TRUE )
                    dtn[dtn$time == s1 & dtn$space == a1 & dtn$technical == g1,'value'] <- sum(true$caa$CAA.disc*cfn, na.rm = TRUE)
                }

                }
            }
        }

        true.dt <- trueData(desc = desc, species = setup.args$species,
                    strataDesc = strIni(timeStrata = switch(paste('s',length(seas), sep = ""), s1 = 'year', s4 = 'quarter', s12 = 'year'),
                                        spaceStrata = ifelse(area == 'all', as.character(NA), 'area') , techStrata = ifelse(gear == 'all', as.character(NA), 'foCatEu5')),
                    laa = laa, lal = lal, ltw =ltw, ltn = ltn, daa = daa,
                    dal = dal, dtw = dtw, dtn = dtn)
      return(true.dt)
}