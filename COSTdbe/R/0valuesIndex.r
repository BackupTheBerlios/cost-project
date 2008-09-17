

library(COSTcore)



#====================================================================
# sampledFO : from a cs consolidated object, a species and a catch category, a vector of the same length as hh columns is returned, with 3 possible values :
#             1 if the haul is sampled/measured and the species is caught in the fraction,
#             0 if the haul is sampled/measured and the species is not caught in the fraction,
#             NA if the haul is not sampled/measured
#====================================================================



extCatchCat <- function(x) {
sapply(x, function(x) substring(as.character(x),1,3))
}


setGeneric("sampledFO", function(x,species,fraction,sampPar=TRUE,...){     #'sampPar' checks if given species is considered as automatically sampled (if TRUE, sppReg=Par <-> sppReg=All, i.e. includes sppReg="Par" in the analysis ) 
	standardGeneric("sampledFO")
	}
)


setMethod("sampledFO", signature(x="csDataCons"), function(x,species,fraction="LAN",sampPar=TRUE,...){

if (missing(species)) stop("Missing species parameter!!")

#Values in 'catReg' field are "All", "Lan", "Dis" & "Non", so...
if (fraction=="DIS") frac <- "Dis" else frac <- "Lan"

  #-----------------------------------------------------------------------------
  # Sampled weight index
  #-----------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    # A haul is considered as sampled (weights) for a given species and a given fraction if :
    #   1) (catchReg=="All") OR (catchReg==frac)
    #   AND
    #   2) (sppReg=="All") OR (sppReg=="Par" AND sampPar=TRUE)
    #---------------------------------------------------------------------------

#hh-index for sampled(1/0)/non sampled(NA) haul (weights) will be the combination of 2 indexes
indexCat <- indexSpp <- rep(0,nrow(x@hh))
#indexCat==1 if catReg=="All" or frac
indexCat[x@hh$catReg%in%c("All",frac)] <- 1
#indexSpp==1 if sppReg=="All" or if sppReg=="Par" & sampPar==TRUE
restrSL <- x@sl[x@sl$spp%in%species & extCatchCat(x@sl$sort)%in%fraction,1:3]
restrSL$ind <- 1 ; indSpp <- merge(x@hh,unique(restrSL[,c(1:2,4)]),all.x=TRUE)$ind     #indSpp <-> index of hauls with related information in sl for given species and fraction
indexSpp[x@hh$sppReg=="All" | (x@hh$sppReg=="Par" & sampPar)] <- 1
#so, Windex = indexCat*indexSpp (sampled haul index)
Windex <- indexCat*indexSpp
indZero <- (Windex==1) & (is.na(indSpp))    #indZero : index of sampled hauls with 0 values 
#finally,...
Windex[Windex==0] <- NA ; Windex[indZero] <- 0


  #-----------------------------------------------------------------------------
  # Sampled numbers index
  #-----------------------------------------------------------------------------

#hh-index for sampled(1/0)/non sampled(NA) hauls (numbers) is Windex, with hauls in sl but not in hl as NA (this means that species was caught in the fraction, but not measured)
#                                                                                                          (O values for weights remain 0 values for numbers)
restrHL <- unique(x@hl[x@hl$spp%in%species & extCatchCat(x@hl$sort)%in%fraction,1:3])
#detect FOs that are recorded in SL but not in HL
restrHL$Ind <- 1 ; indMeas <- merge(unique(restrSL),restrHL,all.x=TRUE) ; indMeas$Ind[is.na(indMeas$Ind)] <- 0
#match index with HH
indMs <- merge(x@hh,indMeas[indMeas$Ind==0,c("PSUid","SSUid","TSUid","Ind")],all.x=TRUE)$Ind
#NAs in 'indMS' means that if info is in SL, then it is in HL
#so, Lindex is...
Lindex <- Windex ; Lindex[!is.na(indMs)] <- NA

#result is returned as a list
return(list(sampWt=Windex,sampLg=Lindex))

})

                  

#--------------------
# Examples
#--------------------

#x <- csDataCons(csDataVal(sole.cs))
#obj <- sampledFO(x,species="Solea vulgaris",fraction="LAN")
#
