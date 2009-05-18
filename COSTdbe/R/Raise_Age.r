#data(sole)
#strD <- strIni(timeStrata="quarter",techStrata="commCat")
#csObject <- csDataCons(csDataVal(subset(sole.cs,sampType%in%c("M","V"))),strD)
#clObject <- clDataCons(clDataVal(sole.cl),strD)
#dbeOutput <- dbeObject(species="Solea solea",catchCat="LAN",strataDesc=strD)
#sex <- as.character(NA)
#
#
#spdAgreg <- COSTcore:::spdAgreg
#
################################################################################
################################################################################
################################################################################
################################################################################



Raise_Age <- function(csObject,dbeOutput,type="fixed",sex=as.character(NA)){#type= "fixed" or "prop" or "ages"

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

#number of fish measured in CA (virtual individuals excluded)
CAreal <- ca[ca$fishId>0,]                                                                              #
nMEAS <- spdAgreg(list(value=CAreal$age),BY=list(time=CAreal$time,space=CAreal$space),function(x) sum(!is.na(x)))       #     
dbeOutput@nMeas$age <- nMEAS                                                                                #

#numbers at length
Ldf <- dbeOutput@lenStruc$estim
N <- tapply(Ldf$value,list(length=Ldf$length,time=Ldf$time,space=Ldf$space,technical=Ldf$technical),sum,na.rm=TRUE)

#creating the stratified ALK with levels from N (duplication for technical strata)
ALK <- tapply(ca$age,list(length=factor(ca$lenCls,levels=dimnames(N)[[1]]),age=ca$age,
                          time=factor(ca$time,levels=dimnames(N)[[2]]),space=factor(ca$space,levels=dimnames(N)[[3]])),length)
ALK[is.na(ALK)] <- 0
  # --> duplication
ll <- dimnames(ALK) ; ll[["technical"]] <- dimnames(N)[[4]]
ALK <- array(rep(as.vector(ALK),dim(N)[4]),dim=c(dim(ALK),dim(N)[4]),dimnames=ll)

  if (all(is.na(dbeOutput@lenStruc$estim))) stop("estimates for length structure are missing in 'dbeOutput' object!!")
#  if (all(is.na(dbeOutput@lenVar))) stop("variance for length structure is missing in 'dbeOutput' object!!")

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
if (!all(apply(Pi.hat,2:4,sum,na.rm=TRUE)==1)) warning("some length classes from 'dbeOutput@lenStruc' slot are not in 'ca' table")

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

VarQij <- aperm(Qij*(1-Qij),c(1,3,4,5,2))/as.vector(nj)
V1 <- aperm(VarQij*as.vector(N*N),c(1,5,2,3,4))

if (!all(is.na(dbeOutput@lenVar))) {
  VarDj <- dbeOutput@lenVar
  VarNj <- tapply(VarDj$value,list(length=VarDj$length,time=VarDj$time,space=VarDj$space,technical=VarDj$technical),sum,na.rm=TRUE)
  V2 <- aperm(aperm(Qij*Qij,c(1,3,4,5,2))*as.vector(VarNj),c(1,5,2,3,4))
  V3 <- aperm(VarQij*as.vector(VarNj),c(1,5,2,3,4))
} else {
  V2 <- V3 <- V1}

#Estimates of total numbers at age
  #total numbers
D.hat <- apply(N,2:4,sum,na.rm=TRUE)
  #total numbers at age
D_i <- Pi.hat*rep(D.hat,each=dim(Pi.hat)[1])

##Estimates of variance at age
#VarDj <- dbeOutput@lenVar
#  #Var(sum(D_j)) = sum(Var(D_j))
#VarD <- tapply(VarDj$value,list(time=factor(VarDj$time,levels=dimnames(N)[[2]]),space=factor(VarDj$space,levels=dimnames(N)[[3]]),
#                                technical=factor(VarDj$technical,levels=dimnames(N)[[4]])),sum,na.rm=TRUE)
#V1 <- Pi.hat*Pi.hat*rep(VarD,each=dim(Pi.hat)[1])
#V2 <- VarPi*rep(D.hat*D.hat,each=dim(VarPi)[1])
#V3 <- VarPi*rep(VarD,each=dim(VarPi)[1])
#
VarD_i <- apply(V1+V2+V3,2:5,sum,na.rm=TRUE)


#results are inserted in dbeOutput object
#####################
  #D_i & VarD_i
df.D_i <- cbind(expand.grid(dimnames(D_i)),value=as.vector(D_i))
df.VarD_i <- cbind(expand.grid(dimnames(VarD_i)),value=as.vector(VarD_i))

df.VarD_i <- df.VarD_i[!is.na(df.D_i$val),] ; df.D_i <- df.D_i[!is.na(df.D_i$val),] 
df.VarD_i <- df.VarD_i[df.D_i$val>0,] ; df.D_i <- df.D_i[df.D_i$val>0,]
  
  #D_i
df.D_i <- df.D_i[order(df.D_i$time,df.D_i$space,df.D_i$technical,df.D_i$age),] ; rownames(df.D_i) <- 1:nrow(df.D_i)
dbeOutput@ageStruc$estim <- df.D_i[,names(dbeOutput@ageStruc$estim)]
  
  #VarD_j
df.VarD_i <- df.VarD_i[order(df.VarD_i$time,df.VarD_i$space,df.VarD_i$technical,df.VarD_i$age),] ; rownames(df.VarD_i) <- 1:nrow(df.VarD_i)

if (!all(is.na(dbeOutput@lenVar))) dbeOutput@ageVar <- df.VarD_i[,names(dbeOutput@ageVar)]

return(dbeOutput)

}



###################
# Exported method #
###################



setGeneric("RaiseAge", function(dbeOutput,
                                 csObject,
                                 type="fixed",
                                 sex=as.character(NA),
                                 ...){
	standardGeneric("RaiseAge")}
)



setMethod("RaiseAge", signature(dbeOutput="dbeOutput",csObject="csDataCons"), function(dbeOutput,
                                                                                       csObject,
                                                                                       type="fixed",
                                                                                       sex=as.character(NA),
                                                                                       ...){
                                                                                                           
Raise_Age(csObject,dbeOutput,type=type,sex=sex)

})





