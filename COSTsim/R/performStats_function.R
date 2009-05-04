PS <- function(estSimObj, trueDataObj, desc, nSamples, ...){

    if(!identical(estSimObj@strataDesc, trueDataObj@strataDesc))
        stop('The stratification in dbeOutputSim and trueData objects  must be the same!')

     desc <- ifelse(missing(desc), estSimObj@desc, desc)

     res <- new('PerformStats', desc = desc, species = trueDataObj@species, strataDesc = trueDataObj@strataDesc,
                catchCat = estSimObj@catchCat, methodDesc = estSimObj@methodDesc, nSamples = nSamples)

     if(estSimObj@catchCat == 'LAN'){  # laa and lal
        res@ageTrue    <- trueDataObj@laa
        res@lenTrue <- trueDataObj@lal

        xx     <- aggregate(list(value = estSimObj@ageStruc$estim$value),
                                    list(time =  estSimObj@ageStruc$estim$time, space = estSimObj@ageStruc$estim$space,
                                         technical = estSimObj@ageStruc$estim$technical, age = estSimObj@ageStruc$estim$age),
                                    mean)
        res@ageEst <- xx[order(as.character(xx$time), as.character(xx$space), as.character(xx$technical), as.numeric(as.character(xx$age))),]
        xx     <- aggregate(list(value = estSimObj@lenStruc$estim$value),
                                    list(time =  estSimObj@lenStruc$estim$time, space = estSimObj@lenStruc$estim$space,
                                         technical = estSimObj@lenStruc$estim$technical, length = estSimObj@lenStruc$estim$length),
                                    mean)
        res@lenEst <- xx[order(as.character(xx$time), as.character(xx$space), as.character(xx$technical), as.numeric(as.character(xx$length))),]

        estj.real.age <- cbind(estSimObj@ageStruc$estim[,-6], value =  estSimObj@ageStruc$estim$value - rep(trueDataObj@laa$value, nSamples))
        estj.real.len <- cbind(estSimObj@lenStruc$estim[,-6], value =  estSimObj@lenStruc$estim$value - rep(trueDataObj@lal$value, nSamples))

        estj.real.age2 <- cbind(estSimObj@ageStruc$estim[,-6], value =  (estSimObj@ageStruc$estim$value - rep(trueDataObj@laa$value, nSamples))^2)
        estj.real.len2 <- cbind(estSimObj@lenStruc$estim[,-6], value =  (estSimObj@lenStruc$estim$value - rep(trueDataObj@lal$value, nSamples))^2)

        estj.real.age.abs <- cbind(estSimObj@ageStruc$estim[,-6], value =  abs(estSimObj@ageStruc$estim$value - rep(trueDataObj@laa$value, nSamples)))
        estj.real.len.abs <- cbind(estSimObj@lenStruc$estim[,-6], value =  abs(estSimObj@lenStruc$estim$value - rep(trueDataObj@lal$value, nSamples)))

        estj.mest.age <- cbind(estSimObj@ageStruc$estim[,-6], value =  estSimObj@ageStruc$estim$value - rep(res@ageEst$value,  nSamples))
        estj.mest.len <- cbind(estSimObj@lenStruc$estim[,-6], value =  estSimObj@lenStruc$estim$value - rep(res@lenEst$value, nSamples))

        estj.real.age.rat <- cbind(estSimObj@ageStruc$estim[,-6], value =  estSimObj@ageStruc$estim$value/rep(trueDataObj@laa$value, nSamples))
        estj.real.len.rat <- cbind(estSimObj@lenStruc$estim[,-6], value =  estSimObj@lenStruc$estim$value/rep(trueDataObj@lal$value, nSamples))

# AGE PERFROMANCE
        sum.estj.real.age <- aggregate(list(value = estj.real.age[,6]), list(time =  estj.real.age$time, space = estj.real.age$space,
                                 technical = estj.real.age$technical, age = estj.real.age$age), sum)
        sum.estj.real.age2 <- aggregate(list(value = estj.real.age2[,6]), list(time =  estj.real.age2$time, space = estj.real.age2$space,
                                 technical = estj.real.age2$technical, age = estj.real.age2$age), sum)
        sum.estj.real.age.abs <- aggregate(list(value = estj.real.age.abs[,6]), list(time =  estj.real.age.abs$time, space = estj.real.age.abs$space,
                                 technical = estj.real.age.abs$technical, age = estj.real.age.abs$age), sum)
        var.estj.real.age <- aggregate(list(value = estj.real.age[,6]), list(time =  estj.real.age$time, space = estj.real.age$space,
                                 technical = estj.real.age$technical, age = estj.real.age$age), var)
        pos.estj.real.age <- aggregate(list(value = ifelse(estj.real.age[,6]>=0, 1, 0)), list(time =  estj.real.age$time, space = estj.real.age$space,
                                 technical = estj.real.age$technical, age = estj.real.age$age), sum)
        sum.estj.real.age.rat <- aggregate(list(value = estj.real.age.rat[,6]*100), list(time =  estj.real.age.rat$time, space = estj.real.age.rat$space,
                                 technical = estj.real.age.rat$technical, age = estj.real.age.rat$age), sum)
        # Bias
        res@ageBias$me  <- cbind(sum.estj.real.age[,-5], value = sum.estj.real.age[,5]/nSamples)
        res@ageBias$poe <- cbind(pos.estj.real.age[,-5], value = pos.estj.real.age[,5]/nSamples)
        res@ageBias$sme <- cbind(sum.estj.real.age[,-5], value = sum.estj.real.age[,5]/(nSamples*trueDataObj@laa$value))
        res@ageBias$par <- cbind(sum.estj.real.age.rat[,-5], value = sum.estj.real.age.rat[,5]/nSamples)
        # Precission
        res@agePrec$var <- var.estj.real.age
        res@agePrec$cv  <- cbind(var.estj.real.age[,-5], value = sqrt(var.estj.real.age[,5])/res@ageEst[,5])
        # Accuracy.
        res@ageAcc$mse <- cbind(sum.estj.real.age2[,-5], value = sum.estj.real.age2[,5]/nSamples)
        res@ageAcc$mae <- cbind(sum.estj.real.age.abs[,-5], value = sum.estj.real.age.abs[,5]/nSamples)
        res@ageAcc$smse <- cbind(sum.estj.real.age2[,-5], value = sum.estj.real.age2[,5]/(nSamples*trueDataObj@laa$value^2))
        res@ageAcc$smae <- cbind(sum.estj.real.age.abs[,-5], value = sum.estj.real.age.abs[,5]/(nSamples*trueDataObj@laa$value))


# LENGTH PERFORMANCE.
        sum.estj.real.len <- aggregate(list(value = estj.real.len[,6]), list(time =  estj.real.len$time, space = estj.real.len$space,
                                 technical = estj.real.len$technical, length = estj.real.len$length), sum)
        sum.estj.real.len2 <- aggregate(list(value = estj.real.len2[,6]), list(time =  estj.real.len2$time, space = estj.real.len2$space,
                                 technical = estj.real.len2$technical, length = estj.real.len2$length), sum)
        sum.estj.real.len.abs <- aggregate(list(value = estj.real.len.abs[,6]), list(time =  estj.real.len.abs$time, space = estj.real.len.abs$space,
                                 technical = estj.real.len.abs$technical, length = estj.real.len.abs$length), sum)
        var.estj.real.len <- aggregate(list(value = estj.real.len[,6]), list(time =  estj.real.len$time, space = estj.real.len$space,
                                 technical = estj.real.len$technical, length = estj.real.len$length), var)
        pos.estj.real.len <- aggregate(list(value = ifelse(estj.real.len[,6]>=0, 1, 0)), list(time =  estj.real.len$time, space = estj.real.len$space,
                                 technical = estj.real.len$technical, length = estj.real.len$length), sum)
        sum.estj.real.len.rat <- aggregate(list(value = estj.real.len.rat[,6]*100), list(time =  estj.real.len.rat$time, space = estj.real.len.rat$space,
                                 technical = estj.real.len.rat$technical, length = estj.real.len.rat$length), sum)
        # Bias
        res@lenBias$me  <- cbind(sum.estj.real.len[,-5], value = sum.estj.real.len[,5]/nSamples)
        res@lenBias$poe <- cbind(pos.estj.real.len[,-5], value = pos.estj.real.len[,5]/nSamples)
        res@lenBias$sme <- cbind(sum.estj.real.len[,-5], value = sum.estj.real.len[,5]/(nSamples*trueDataObj@lal$value))
        res@lenBias$par <- cbind(sum.estj.real.len.rat[,-5], value = sum.estj.real.len.rat[,5]/nSamples)
        # Precission
        res@lenPrec$var <- var.estj.real.len
        res@lenPrec$cv  <- cbind(var.estj.real.len[,-5], value = sqrt(var.estj.real.len[,5])/res@lenEst[,5])
        # Accuracy.
        res@lenAcc$mse <- cbind(sum.estj.real.len2[,-5], value = sum.estj.real.len2[,5]/nSamples)
        res@lenAcc$mae <- cbind(sum.estj.real.len.abs[,-5], value = sum.estj.real.len.abs[,5]/nSamples)
        res@lenAcc$smse <- cbind(sum.estj.real.len2[,-5], value = sum.estj.real.len2[,5]/(nSamples*trueDataObj@lal$value^2))
        res@lenAcc$smae <- cbind(sum.estj.real.len.abs[,-5], value = sum.estj.real.len.abs[,5]/(nSamples*trueDataObj@lal$value))
    }
    else{ # DIS
        ageTrue    <- trueDataObj@laa
        lenTrue <- trueDataObj@lal
    }
  return(res)
  }