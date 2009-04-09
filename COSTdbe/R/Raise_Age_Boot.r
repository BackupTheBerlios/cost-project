# Raise_Age_Boot

## FUNCTION STILL IN DEVELOPMENT,
#need to check iters that don't give a particular age-length combination are assigned 0 for calculating mean and var across iter

#library(COSTcore) #1.3-3
#library(COSTdbe)  #0.1-9

#spdAgreg <- COSTcore:::spdAgreg
#resample <- COSTdbe:::resample
#As.num <- COSTcore:::As.num

#data(sole)
#strD <- strIni(timeStrata="quarter",techStrata="commCat")
#csObject <- csDataCons(csDataVal(subset(sole.cs,sampType%in%c("M","V"))),strD)
#clObject <- clDataCons(clDataVal(sole.cl),strD)
#dbeOutput <- dbeObject(species="Solea solea",catchCat="LAN",strataDesc=strD)
#sex <- as.character(NA)
#
#sol.dbe.an <- RaiseLgth (dbeOutput = dbeOutput, csObject = csObject, clObject = clObject)
#sol.dbe.boot <- RaiseLgthBoot (dbeOutput = dbeOutput, csObject = csObject, clObject = clObject, B=10)

#sol.dbe.an <- RaiseAge (csObject = csObject, dbeOutput = sol.dbe.an, type="fixed")
#sol.dbe.boot <- RaiseAgeBoot (csObject = csObject, dbeOutput = sol.dbe.boot, type="fixed")


################################################################################
################################################################################
################################################################################
################################################################################



Raise_Age_Boot <- function(csObject,dbeOutput,type="fixed",sex=as.character(NA)){#type= "fixed" or "prop" or "ages"

sp <- dbeOutput@species
ca <- ca(csObject)
#ca table is subset
ca <- ca[ca$spp%in%sp,]   
if (nrow(ca)==0) stop("no CA data for specified species in input object!!")

#number of age samples (PSUid/SSUid) is calculated at this stage (before subsetting on 'sex')               #
Unit <- paste(ca$PSUid,ca$SSUid,sep=":-:")                                                                  #
nSAMP <- spdAgreg(list(value=Unit),BY=list(time=ca$time,space=ca$space),function(x) length(unique(x)))      #
dbeOutput@nSamp$age <- nSAMP                                                                                #  ADDED : MM 02/04/2009

                                                                                          #

if (!(all(is.na(sex)))) {ca <- ca[ca$sex%in%sex,]                                                            
                         if (nrow(ca)==0) stop("no CA data for specified sex in input object!!")            #                                                                               #
}                                                                                                           #

#number of fish measured in HL                                                                              #
nMEAS <- spdAgreg(list(value=ca$age),BY=list(time=ca$time,space=ca$space),function(x) sum(!is.na(x)))       #     
dbeOutput@nMeas$age <- nMEAS                                                                                #

#numbers at length all replicates
if (all(is.na(dbeOutput@lenStruc$rep))) stop("estimates for length structure are missing in 'dbeOutput' object!!")

Ldfall <- dbeOutput@lenStruc$rep
# set number of iterations to match reps of length
B <- max(Ldfall$iter)

# PSUids for each strata combination
CASTR <- paste(ca$time,ca$space,sep=":-:")
UnitSTR = data.frame (Unit, CASTR, stringsAsFactors = F)

ageids <- unique( UnitSTR )
ageids <- ageids[order(ageids$CASTR),] # ageids and nSAMP should now have STR in same order
uStrat <- paste (nSAMP$time,nSAMP$space,sep=":-:")
#identify start and end positions for each strata
start.pos = c(1, (cumsum(nSAMP$value)[-length(nSAMP$value)]) +1 )
end.pos = cumsum(nSAMP$value)

# sample new set of PSUid for each strata combination - stratified bootstrap resampling - may be able to use boot function instead
bootAgeid = matrix (NA, nrow = dim(ageids)[1], ncol=B+1)
dimnames(bootAgeid) = list(NULL, c("orig", paste("iter.",1:B,sep="")))
# original sample ids in first column
bootAgeid[,1] = ageids$Unit

# put resampled ids for each strata into vector in relevant places,
# assigning to all columns using size = nSAMP$value[i] * B, instead of using loop by iteration 1 to B
for ( i in 1:length(uStrat) ){
# order of nSAMP$value and uStrat needs to match
  bootAgeid [ start.pos[i]:end.pos[i], -1] = resample (ageids$Unit [ ageids$CASTR == uStrat[i] ], size = nSAMP$value[i] * B, replace=T )
  }

CA.orig <- ca
CA.orig$Unit <- paste(CA.orig$PSUid,CA.orig$SSUid,sep=":-:")


####### START OF BOOTSTRAP LOOP #################
ac.list = vector("list", B+1)

# i=1 uses original data, its output is labelled iter=0
for (i in 1:(B+1) ){
print(i-1)

ca = merge(data.frame(Unit = bootAgeid[,i]), CA.orig, by="Unit")

#numbers at length
Ldf <- Ldfall[Ldfall$iter == (i-1),]
N <- tapply(Ldf$value,list(length=Ldf$length,time=Ldf$time,space=Ldf$space,technical=Ldf$technical),sum,na.rm=TRUE)

#creating the stratified ALK with levels from N (duplication for technical strata)
ALK <- tapply(ca$age,list(length=factor(ca$lenCls,levels=dimnames(N)[[1]]),age=ca$age,
                          time=factor(ca$time,levels=dimnames(N)[[2]]),space=factor(ca$space,levels=dimnames(N)[[3]])),length)
ALK[is.na(ALK)] <- 0

## THIS SECTION IS COPIED FROM VesselRaiseBoot as a reminder to do something about gaps in alks
# to use it will need to setup alkgaps.list and alkgaps.counter and change alk to ALK
## Fill in missing rows of alk  (Will be quicker to only run this if there are gaps that need filling)
# (all gaps of length <= value are filled with the sum of surrounding filled classes)

# need to have functions in alkGaps.r available
# Disadvantage is that it creates 'virtual' otoliths

#gaps.out = gapsRm(alk,type="fillMiss",value=2,preview=FALSE,postview=FALSE)
#alk = gaps.out$alk
# records with otoliths added
#if ( dim(gaps.out$addIndTab)[1]>0 ){
#alkgaps.list[[i]] = cbind(gaps.out$addIndTab, iter=(i-1))
#alkgaps.counter = alkgaps.counter + 1 }


  # --> duplication
ll <- dimnames(ALK) ; ll[["technical"]] <- dimnames(N)[[4]]
ALK <- array(rep(as.vector(ALK),dim(N)[4]),dim=c(dim(ALK),dim(N)[4]),dimnames=ll)

#nj : number of sampled individuals of length j
  nj <- tapply(ca$lenCls,list(length=factor(ca$lenCls,levels=dimnames(N)[[1]]),time=factor(ca$time,levels=dimnames(N)[[2]]),
                              space=factor(ca$space,levels=dimnames(N)[[3]])),function(x) sum(!is.na(x)))
  nj[is.na(nj)] <- 0
  # --> duplication
  ll2 <- dimnames(nj) ; ll2[["technical"]] <- dimnames(N)[[4]]
  nj <- array(rep(as.vector(nj),dim(N)[4]),dim=c(dim(nj),dim(N)[4]),dimnames=ll2)

#Nl : number of length-sampled individuals
	Nl <- apply(nj,2:4,sum)
#ns : number of samples to be aged
	ns <- apply(ALK,3:5,sum)
#nl : number of sampled length categories
  nl <- apply(nj,2:4,function(x) sum(x>0))
#Q'ij
	Qij <- aperm(aperm(ALK,c(1,3,4,5,2))/as.vector(apply(ALK,c(1,3:5),sum)),c(1,5,2:4))
#njStar
	if (type=="fixed") {
    njStar <- ns/nl
  } else {
    njStar <- nj*rep(as.vector(ns/Nl),each=dim(nj)[1])}
#lj : proportion of length j in the population
	lj <- N/rep(as.vector(apply(N,2:4,sum,na.rm=TRUE)),each=dim(N)[1])
  lj[is.na(lj)] <- 0

#'pi' calculation
#------------
Pi.hat <- apply(aperm(aperm(Qij,c(1,3,4,5,2))*as.vector(lj),c(1,5,2,3,4)),2:5,sum,na.rm=TRUE)
#test on Pi.hat (missing length class in ca)
#if (!all(apply(Pi.hat,2:4,sum,na.rm=TRUE)==1)) warning("it seems that some length classes from 'dbeOutput@lenStruc' slot are not in 'ca' table")

#Var.pi calculation
#----------------
#if (type=="ages") {
#	a1 <- Pi.hat*(1-Pi.hat)
#	VarPi <- a1/rep(ns,each=dim(a1)[1])
#} else {
#if (type=="fixed") {
#	b1 <- apply(aperm(aperm(Qij*(1-Qij),c(1,3,4,5,2))*as.vector(lj*(1-lj)),c(1,5,2,3,4)),2:5,sum,na.rm=TRUE)
#	b2 <- apply(aperm(aperm(Qij*(1-Qij),c(1,3,4,5,2))*as.vector(lj*lj),c(1,5,2,3,4)),2:5,sum,na.rm=TRUE)
#	b3 <- apply(aperm(aperm(Qij*Qij,c(1,3,4,5,2))*as.vector(lj),c(1,5,2,3,4)),2:5,sum,na.rm=TRUE)-Pi.hat^2
#  VarPi <- b1/rep(Nl*njStar,each=dim(b1)[1]) + b2/rep(njStar,each=dim(b2)[1]) + b3/rep(Nl,each=dim(b3)[1])
#} else {   #i.e if (type=="prop")
#	c1 <- apply(aperm(aperm(Qij*(1-Qij),c(1,3,4,5,2))*as.vector(lj),c(1,5,2,3,4)),2:5,sum,na.rm=TRUE)
#	c2 <- apply(aperm(aperm(Qij*Qij,c(1,3,4,5,2))*as.vector(lj),c(1,5,2,3,4)),2:5,sum,na.rm=TRUE)-Pi.hat^2
#  VarPi <- c1/rep(ns,each=dim(c1)[1]) + c2/rep(Nl,each=dim(c2)[1])
#}}

#Estimates of total numbers at age
  #total numbers
D.hat <- apply(N,2:4,sum,na.rm=TRUE)
  #total numbers at age
D_i <- Pi.hat*rep(D.hat,each=dim(Pi.hat)[1])

#Estimates of variance at age
#VarDj <- dbeOutput@lenVar
  #Var(sum(D_j)) = sum(Var(D_j))
#VarD <- tapply(VarDj$value,list(time=factor(VarDj$time,levels=dimnames(N)[[2]]),space=factor(VarDj$space,levels=dimnames(N)[[3]]),
#                                technical=factor(VarDj$technical,levels=dimnames(N)[[4]])),sum,na.rm=TRUE)
#V1 <- Pi.hat*Pi.hat*rep(VarD,each=dim(Pi.hat)[1])
#V2 <- VarPi*rep(D.hat*D.hat,each=dim(VarPi)[1])
#V3 <- VarPi*rep(VarD,each=dim(VarPi)[1])

#VarD_i <- V1+V2+V3


#results are inserted in dbeOutput object
#####################
  #D_i & VarD_i
df.D_i <- cbind(expand.grid(dimnames(D_i)),value=as.vector(D_i))
#df.VarD_i <- cbind(expand.grid(dimnames(VarD_i)),value=as.vector(VarD_i))

#df.VarD_i <- df.VarD_i[!is.na(df.D_i$val),] ; df.D_i <- df.D_i[!is.na(df.D_i$val),]
#df.VarD_i <- df.VarD_i[df.D_i$val>0,] ; df.D_i <- df.D_i[df.D_i$val>0,]
  
  #D_i
df.D_i <- df.D_i[order(df.D_i$time,df.D_i$space,df.D_i$technical,df.D_i$age),] ; rownames(df.D_i) <- 1:nrow(df.D_i)
ac.list[[i]] <- df.D_i[,names(dbeOutput@ageStruc$estim)] # changed for boot

  #VarD_j
#df.VarD_i <- df.VarD_i[order(df.VarD_i$time,df.VarD_i$space,df.VarD_i$technical,df.VarD_i$age),] ; rownames(df.VarD_i) <- 1:nrow(df.VarD_i)
#dbeOutput@ageVar <- df.VarD_i[,names(dbeOutput@ageVar)]

} # End of bootstrap loop

# Age structure
#convert list of length distributions into data.frame matching dbeOutput format
rep.iter = unlist (lapply(ac.list, FUN = function(x) {dim(x)[1]}) )

ac.df = dbeOutput@ageStruc$rep = data.frame  (time =  unlist(lapply(ac.list, FUN = function(x) {x[,"time"]})),
                                space =  unlist(lapply(ac.list, FUN = function(x) {x[,"space"]})),
                                technical = unlist(lapply(ac.list, FUN = function(x) {x[,"technical"]})),
                                age =  unlist(lapply(ac.list, FUN = function(x) {x[,"age"]})),
                                value =  unlist(lapply(ac.list, FUN = function(x) {x[,"value"]})),
                                iter =  rep(0:B, times = rep.iter)   )

ac.df = ac.df[ac.df$iter > 0,]
ac.mean = spdAgreg (list (value=ac.df$value), BY = list(time=ac.df$time, space=ac.df$space, technical=ac.df$technical, age=ac.df$age), mean)
ac.mean$age = As.num(ac.mean$age)
ac.mean = ac.mean [order(ac.mean$time, ac.mean$space, ac.mean$technical, ac.mean$age),]
dimnames(ac.mean)[[1]] = 1:(dim(ac.mean)[1])

ac.var = spdAgreg (list (value=ac.df$value), BY = list(time=ac.df$time, space=ac.df$space, technical=ac.df$technical, age=ac.df$age), var)
ac.var$age = As.num(ac.var$age)
ac.var = ac.var [order(ac.var$time, ac.var$space, ac.var$technical, ac.var$age),]
dimnames(ac.var)[[1]] = 1:(dim(ac.var)[1])

dbeOutput@ageStruc$estim = ac.mean
dbeOutput@ageVar = ac.var

return(dbeOutput)

}



###################
# Exported method #
###################



setGeneric("RaiseAgeBoot", function(dbeOutput,
                                 csObject,
                                 type="fixed",
                                 sex=as.character(NA),
                                 ...){
	standardGeneric("RaiseAgeBoot")}
)


setMethod("RaiseAgeBoot", signature(dbeOutput="dbeOutput",csObject="csDataCons"), function(dbeOutput,
                                                                                       csObject,
                                                                                       type="fixed",
                                                                                       sex=as.character(NA),
                                                                                       ...){
                                                                                                           
Raise_Age_Boot(csObject=csObject,dbeOutput=dbeOutput,type=type,sex=sex)

})





