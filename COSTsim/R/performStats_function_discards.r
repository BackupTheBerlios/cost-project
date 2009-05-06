PS.dis <- function(estSimObj, trueDataObj, desc, nSamples){

        res <- new('PerformStats', desc = desc, species = trueDataObj@species, strataDesc = trueDataObj@strataDesc,
                catchCat = estSimObj@catchCat, methodDesc = estSimObj@methodDesc, nSamples = nSamples)

        res@totalWTrue     <- trueDataObj@dtw
        res@lenTrue        <- trueDataObj@dal

        xx     <- aggregate(list(value = estSimObj@totalW$estim$value),
                                    list(time =  estSimObj@totalW$estim$time, space = estSimObj@totalW$estim$space,
                                         technical = estSimObj@totalW$estim$technical),
                                    mean)
        res@totalWEst <- xx[order(as.character(xx$time), as.character(xx$space), as.character(xx$technical)),]
        xx     <- aggregate(list(value = estSimObj@lenStruc$estim$value),
                                    list(time =  estSimObj@lenStruc$estim$time, space = estSimObj@lenStruc$estim$space,
                                         technical = estSimObj@lenStruc$estim$technical, length = estSimObj@lenStruc$estim$length),
                                    mean)
        res@lenEst <- xx[order(as.character(xx$time), as.character(xx$space), as.character(xx$technical), as.numeric(as.character(xx$length))),]

        estj.real.w   <- cbind(estSimObj@totalW$estim[,-5], value =  estSimObj@totalW$estim$value - rep(trueDataObj@dtw$value, nSamples))
        estj.real.len <- cbind(estSimObj@lenStruc$estim[,-6], value =  estSimObj@lenStruc$estim$value - rep(trueDataObj@dal$value, nSamples))

        estj.real.w2   <- cbind(estSimObj@totalW$estim[,-5], value =  (estSimObj@totalW$estim$value - rep(trueDataObj@dtw$value, nSamples))^2)
        estj.real.len2 <- cbind(estSimObj@lenStruc$estim[,-6], value =  (estSimObj@lenStruc$estim$value - rep(trueDataObj@dal$value, nSamples))^2)

        estj.real.w.abs   <- cbind(estSimObj@totalW$estim[,-5], value =  abs(estSimObj@totalW$estim$value - rep(trueDataObj@dtw$value, nSamples)))
        estj.real.len.abs <- cbind(estSimObj@lenStruc$estim[,-6], value =  abs(estSimObj@lenStruc$estim$value - rep(trueDataObj@dal$value, nSamples)))

        estj.mest.w   <- cbind(estSimObj@totalW$estim[,-5], value =  estSimObj@totalW$estim$value - rep(res@totalWEst$value,  nSamples))
        estj.mest.len <- cbind(estSimObj@lenStruc$estim[,-6], value =  estSimObj@lenStruc$estim$value - rep(res@lenEst$value, nSamples))

        estj.real.w.rat   <- cbind(estSimObj@totalW$estim[,-5], value =  estSimObj@totalW$estim$value/rep(trueDataObj@dtw$value, nSamples))
        estj.real.len.rat <- cbind(estSimObj@lenStruc$estim[,-6], value =  estSimObj@lenStruc$estim$value/rep(trueDataObj@dal$value, nSamples))

# TotalW PERFROMANCE
        sum.estj.real.w       <- aggregate(list(value = estj.real.w[,5]), list(time =  estj.real.w$time, space = estj.real.w$space,
                                        technical = estj.real.w$technical), sum)
        sum.estj.real.w2      <- aggregate(list(value = estj.real.w2[,5]), list(time =  estj.real.w2$time, space = estj.real.w2$space,
                                        technical = estj.real.w2$technical), sum)
        sum.estj.real.w.abs   <- aggregate(list(value = estj.real.w.abs[,5]), list(time =  estj.real.w.abs$time, space = estj.real.w.abs$space,
                                        technical = estj.real.w.abs$technical), sum)
        var.estj.real.w       <- aggregate(list(value = estj.real.w[,5]), list(time =  estj.real.w$time, space = estj.real.w$space,
                                        technical = estj.real.w$technical), var)
        pos.estj.real.w       <- aggregate(list(value = ifelse(estj.real.w[,5]>=0, 1, 0)), list(time =  estj.real.w$time, space = estj.real.w$space,
                                        technical = estj.real.w$technical), sum)
        sum.estj.real.w.rat   <- aggregate(list(value = estj.real.w.rat[,5]*100), list(time =  estj.real.w.rat$time, space = estj.real.w.rat$space,
                                        technical = estj.real.w.rat$technical), sum)
       
        sum.estj.real.w       <- order.df.disc(sum.estj.real.w, len = FALSE)
        sum.estj.real.w2      <- order.df.disc(sum.estj.real.w2, len = FALSE)
        sum.estj.real.w.abs   <- order.df.disc(sum.estj.real.w.abs, len = FALSE)
        var.estj.real.w       <- order.df.disc(var.estj.real.w, len = FALSE)
        pos.estj.real.w       <- order.df.disc(pos.estj.real.w, len = FALSE)
        sum.estj.real.w.rat   <- order.df.disc(sum.estj.real.w.rat, len = FALSE)
        
        
        # Bias
        res@totalWBias$me  <- cbind(sum.estj.real.w[,-4], value = sum.estj.real.w[,4]/nSamples)
        res@totalWBias$poe <- cbind(pos.estj.real.w[,-4], value = pos.estj.real.w[,4]/nSamples)
        res@totalWBias$sme <- cbind(sum.estj.real.w[,-4], value = sum.estj.real.w[,4]/(nSamples*trueDataObj@dtw$value))
        res@totalWBias$par <- cbind(sum.estj.real.w.rat[,-4], value = sum.estj.real.w.rat[,4]/nSamples)
        # Precission
        res@totalWPrec$var <- var.estj.real.w
        res@totalWPrec$cv  <- cbind(var.estj.real.w[,-4], value = sqrt(var.estj.real.w[,4])/res@totalWEst[,4])
        # Accuracy.
        res@totalWAcc$mse  <- cbind(sum.estj.real.w2[,-4], value = sum.estj.real.w2[,4]/nSamples)
        res@totalWAcc$mae  <- cbind(sum.estj.real.w.abs[,-4], value = sum.estj.real.w.abs[,4]/nSamples)
        res@totalWAcc$smse <- cbind(sum.estj.real.w2[,-4], value = sum.estj.real.w2[,4]/(nSamples*trueDataObj@dtw$value^2))
        res@totalWAcc$smae <- cbind(sum.estj.real.w.abs[,-4], value = sum.estj.real.w.abs[,4]/(nSamples*trueDataObj@dtw$value))


# LENGTH PERFORMANCE.
        sum.estj.real.len       <- aggregate(list(value = estj.real.len[,6]), list(time =  estj.real.len$time, space = estj.real.len$space,
                                        technical = estj.real.len$technical, length = estj.real.len$length), sum)
        sum.estj.real.len2      <- aggregate(list(value = estj.real.len2[,6]), list(time =  estj.real.len2$time, space = estj.real.len2$space,
                                        technical = estj.real.len2$technical, length = estj.real.len2$length), sum)
        sum.estj.real.len.abs   <- aggregate(list(value = estj.real.len.abs[,6]), list(time =  estj.real.len.abs$time, space = estj.real.len.abs$space,
                                        technical = estj.real.len.abs$technical, length = estj.real.len.abs$length), sum)
        var.estj.real.len       <- aggregate(list(value = estj.real.len[,6]), list(time =  estj.real.len$time, space = estj.real.len$space,
                                        technical = estj.real.len$technical, length = estj.real.len$length), var)
        pos.estj.real.len       <- aggregate(list(value = ifelse(estj.real.len[,6]>=0, 1, 0)), list(time =  estj.real.len$time, space = estj.real.len$space,
                                        technical = estj.real.len$technical, length = estj.real.len$length), sum)
        sum.estj.real.len.rat   <- aggregate(list(value = estj.real.len.rat[,6]*100), list(time =  estj.real.len.rat$time, space = estj.real.len.rat$space,
                                        technical = estj.real.len.rat$technical, length = estj.real.len.rat$length), sum)

        sum.estj.real.len       <- order.df.disc(sum.estj.real.len, len = TRUE)
        sum.estj.real.len2      <- order.df.disc(sum.estj.real.len2, len = TRUE)
        sum.estj.real.len.abs   <- order.df.disc(sum.estj.real.len.abs, len = TRUE)
        var.estj.real.len       <- order.df.disc(var.estj.real.len, len = TRUE)
        pos.estj.real.len       <- order.df.disc(pos.estj.real.len, len = TRUE)
        sum.estj.real.len.rat   <- order.df.disc(sum.estj.real.len.rat, len = TRUE)
        
        # Bias
        res@lenBias$me  <- cbind(sum.estj.real.len[,-5], value = sum.estj.real.len[,5]/nSamples)
        res@lenBias$poe <- cbind(pos.estj.real.len[,-5], value = pos.estj.real.len[,5]/nSamples)
        res@lenBias$sme <- cbind(sum.estj.real.len[,-5], value = sum.estj.real.len[,5]/(nSamples*trueDataObj@dal$value))
        res@lenBias$par <- cbind(sum.estj.real.len.rat[,-5], value = sum.estj.real.len.rat[,5]/nSamples)
        # Precission
        res@lenPrec$var <- var.estj.real.len
        res@lenPrec$cv  <- cbind(var.estj.real.len[,-5], value = sqrt(var.estj.real.len[,5])/res@lenEst[,5])
        # Accuracy.
        res@lenAcc$mse  <- cbind(sum.estj.real.len2[,-5], value = sum.estj.real.len2[,5]/nSamples)
        res@lenAcc$mae  <- cbind(sum.estj.real.len.abs[,-5], value = sum.estj.real.len.abs[,5]/nSamples)
        res@lenAcc$smse <- cbind(sum.estj.real.len2[,-5], value = sum.estj.real.len2[,5]/(nSamples*trueDataObj@dal$value^2))
        res@lenAcc$smae <- cbind(sum.estj.real.len.abs[,-5], value = sum.estj.real.len.abs[,5]/(nSamples*trueDataObj@dal$value))

  return(res)
}
  
order.df.disc <- function(df, len = TRUE){
     if(len == TRUE)
         df <- df[order(as.character(df$time), as.character(df$space), as.character(df$technical), as.numeric(as.character(df$length))),]
     else
        df <- df[order(as.character(df$time), as.character(df$space), as.character(df$technical)),]

        return(df)
}