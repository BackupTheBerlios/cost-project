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



Raise_Lgth <- function(csObject,clObject,dbeOutput,spp,taxon,sex=as.character(NA)){

sp <- dbeOutput@species
if (missing(taxon)) taxon <- sp
if (missing(spp)) spp <- sp   # todo : give as an argument the species in CS that correspond to taxon in CL

# weight of the species/taxa/sex in the sampled box
# this weight is to be calculated by summing the sampled weight by species/taxa/sex before subseting
# and added as a third weight variable

indNewWt <- TRUE
if (is.na(sex) & all(taxon%in%sp)) indNewWt <- FALSE
if (indNewWt) {

  NW <- spdAgreg(list(newWt=csObject@sl$subSampWt),BY=list(PSUid=csObject@sl$PSUid,SSUid=csObject@sl$SSUid,TSUid=csObject@sl$TSUid,
                                                         time=csObject@sl$time,space=csObject@sl$space,
                                                         technical=csObject@sl$technical,sort=csObject@sl$sort),sum,na.rm=TRUE)

  mergeSL <- merge(csObject@sl,NW,all.x=TRUE,sort=FALSE)
  ssw <- csObject@sl$subSampWt
  csObject@sl$subSampWt <- mergeSL$newWt
  csObject@sl$lenCode <- ssw
} else {
  csObject@sl$lenCode <- csObject@sl$subSampWt
}


#subsetting "csObject" 
csObject <- subsetSpp(csObject,spp%in%sp)
if (!is.na(sex)) csObject <- subsetSpp(csObject,sex%in%sex)  #to keep all 'tr' and 'hh' information, we use "subsetSpp" method
csObject <- subsetSpp(csObject,spp%in%sp)
#subsetting "clObject"
clObject <- subset(clObject,taxon%in%taxon)    

SL <- csObject@sl
HL <- csObject@hl
#creating PSTUid field, concatenation of PSUid, SSUid, TSUid
SL$PSTUid <- apply(SL[,c("PSUid","SSUid","TSUid")],1,function(x) paste(x,collapse=":-:"))
HL$PSTUid <- apply(HL[,c("PSUid","SSUid","TSUid")],1,function(x) paste(x,collapse=":-:"))
#creating STR field, concatenation of time, space, technical
HHSTR <- apply(csObject@hh[,c("time","space","technical")],1,function(x) paste(x,collapse=":-:"))
SL$STR <- apply(SL[,c("time","space","technical")],1,function(x) paste(x,collapse=":-:"))
HL$STR <- apply(HL[,c("time","space","technical")],1,function(x) paste(x,collapse=":-:"))

if (nrow(SL)==0) stop("the arguments given resulted in an empty table!!")

#total landed weights per strata
clObject@cl$landMult[is.na(clObject@cl$landMult)] <- 1
#TotLand = OffLand*Multi + UnallocCat + MisallocCat
totLand <- mapply(function(w,x,y,z) sum(c(w*x,y,z),na.rm=TRUE),clObject@cl$landWt,clObject@cl$landMult,clObject@cl$unallocCatchWt,clObject@cl$misRepCatchWt)
totLandings <- spdAgreg(list(W=totLand),BY=list(time=clObject@cl$time,space=clObject@cl$space,technical=clObject@cl$technical),sum,na.rm=TRUE)


#weight of the level
wl <- tapply(SL$wt,list(STR=SL$STR,sort=SL$sort,TSUid=SL$TSUid,SSUid=SL$SSUid,PSUid=SL$PSUid),sum,na.rm=TRUE)   #reference for factor levels
#sampled weight
ws <- tapply(SL$subSampWt,list(STR=SL$STR,sort=SL$sort,TSUid=SL$TSUid,SSUid=SL$SSUid,PSUid=SL$PSUid),sum,na.rm=TRUE)
#weight of the species/taxa/sex
wt  <- tapply(as.numeric(as.character(SL$lenCode)),list(STR=SL$STR,sort=SL$sort,TSUid=SL$TSUid,SSUid=SL$SSUid,PSUid=SL$PSUid),sum,na.rm=TRUE)

#number of fish in the sample by length
d_j <- tapply(HL$lenNum,list(STR=factor(HL$STR,levels=levels(factor(SL$STR))),
                              sort=factor(HL$sort,levels=levels(factor(SL$sort))),
                              TSUid=factor(HL$TSUid,levels=levels(factor(SL$TSUid))),
                              SSUid=factor(HL$SSUid,levels=levels(factor(SL$SSUid))),
                              PSUid=factor(HL$PSUid,levels=levels(factor(SL$PSUid))),
                              lenCls=HL$lenCls),sum,na.rm=TRUE)

#TSUid stage
d_jtsu <- apply(d_j*(as.vector(wl/ws)),c(1,3:6),sum,na.rm=TRUE)
w_tsu <- apply(wt*wl/ws,c(1,3:5),sum,na.rm=TRUE)  # A REVOIR
wl_tsu <- apply(wl,c(1,3:5),sum,na.rm=TRUE)

#SSUid stage
sum.d_jssu <- apply(d_jtsu,c(1,4,5),sum,na.rm=TRUE)
sum.w_ssu <- apply(w_tsu,c(1,4),sum,na.rm=TRUE)
sum.wl_ssu <- apply(wl_tsu,c(1,4),sum,na.rm=TRUE)

#number of SSUid (total and sampled)
Mi <- tapply(csObject@hh$SSUid,list(STR=factor(HHSTR,levels=levels(factor(SL$STR))),
                                      PSUid=factor(csObject@hh$PSUid,levels=levels(factor(SL$PSUid)))),function(x) length(unique(x)))
mi <- tapply(HL$SSUid,list(STR=factor(HL$STR,levels=levels(factor(SL$STR))),
                           PSUid=factor(HL$PSUid,levels=levels(factor(SL$PSUid)))),function(x) length(unique(x)))

#PSUid stage                           
d_jpsu <- sum.d_jssu*as.vector(Mi/mi)
w_psu <- sum.w_ssu*Mi/mi
wl_psu <- sum.wl_ssu*Mi/mi

#sum.wt <- apply(wt,1,sum,na.rm=TRUE)
#sum.ws <- apply(ws,1,sum,na.rm=TRUE)

sum.d_jpsu <- apply(d_jpsu,c(1,3),sum,na.rm=TRUE)
sum.w_psu <- apply(w_psu,1,sum,na.rm=TRUE)
sum.wl_psu <- apply(wl_psu,1,sum,na.rm=TRUE)

W <- tapply(totLandings$W*1000,list(factor(apply(totLandings[,c("time","space","technical")],1,function(x) paste(x,collapse=":-:")),levels=dimnames(d_j)[[1]])),sum,na.rm=TRUE)

D_j <- sum.d_jpsu*as.vector(W/sum.wl_psu)
WHat <-  sum.w_psu*W/sum.wl_psu
                                                                                                    


## variance calculation based on 'Detecting sampling outliers and sampling heterogeneity when...' article (J. Vigneau & S. Mah�vas)
######################

#n : number of sampled PSUid per strata
n <- tapply(HL$PSUid,list(STR=factor(HL$STR,levels=levels(factor(SL$STR)))),function(x) length(unique(x)))

first  <- W^2
second <- (1-(sum.wl_psu/W))/(sum.wl_psu^2/n) 
  third.1  <- sum.d_jpsu/as.vector(sum.wl_psu)
  third.2 <- aperm(array(rep(as.vector(third.1),dim(w_psu)[2]),dim=dim(d_jpsu)[c(1,3,2)]),c(1,3,2))
  third.3 <- d_jpsu-third.2*as.vector(wl_psu)
third <- apply(third.3^2,c(1,3),sum,na.rm=TRUE)/as.vector(n-1)
VarD_j <- third*as.vector(first*second)



#results are inserted in dbeOutput object
#####################
  #D_j & VarD_j
df.D_j <- cbind(expand.grid(dimnames(D_j)),value=as.vector(D_j))
df.VarD_j <- cbind(expand.grid(dimnames(VarD_j)),value=as.vector(VarD_j))

df.VarD_j <- df.VarD_j[!is.na(df.D_j$val),] ; df.D_j <- df.D_j[!is.na(df.D_j$val),] 
df.VarD_j <- df.VarD_j[df.D_j$val>0,] ; df.D_j <- df.D_j[df.D_j$val>0,]
  
  #D_j
df.D_j <- cbind(df.D_j,do.call("rbind",lapply(as.character(df.D_j$STR),function(x) strsplit(x,":-:")[[1]])))
names(df.D_j) <- c("STR","length","value","time","space","technical")
df.D_j <- df.D_j[order(df.D_j$time,df.D_j$space,df.D_j$technical,df.D_j$length),] ; rownames(df.D_j) <- 1:nrow(df.D_j)
dbeOutput@lenStruc$estim <- df.D_j[,names(dbeOutput@lenStruc$estim)]
  
  #VarD_j
df.VarD_j <- cbind(df.VarD_j,do.call("rbind",lapply(as.character(df.VarD_j$STR),function(x) strsplit(x,":-:")[[1]])))
names(df.VarD_j) <- c("STR","length","value","time","space","technical")
df.VarD_j <- df.VarD_j[order(df.VarD_j$time,df.VarD_j$space,df.VarD_j$technical,df.VarD_j$length),] ; rownames(df.VarD_j) <- 1:nrow(df.VarD_j)
dbeOutput@lenVar <- df.VarD_j[,names(dbeOutput@lenVar)]

  #WHat
df.WHat <- cbind(value=WHat,as.data.frame(do.call("rbind",lapply(names(WHat),function(x) strsplit(x,":-:")[[1]]))))
names(df.WHat) <- c("value","time","space","technical")
df.WHat <- df.WHat[order(df.WHat$time,df.WHat$space,df.WHat$technical),] ; rownames(df.WHat) <- 1:nrow(df.WHat)
dbeOutput@totalW$estim <- df.WHat[,names(dbeOutput@totalW$estim)]

return(dbeOutput)

}

###################
# Exported method #
###################



setGeneric("RaiseLgth", function(dbeOutput,
                                 csObject,
                                 clObject,
                                 spp,
                                 taxon,
                                 sex=as.character(NA),
                                 ...){
	standardGeneric("RaiseLgth")}
)



setMethod("RaiseLgth", signature(dbeOutput="dbeOutput",csObject="csDataCons",clObject="clDataCons"), function(dbeOutput,
                                                                                                              csObject,
                                                                                                              clObject,
                                                                                                              spp,
                                                                                                              taxon,
                                                                                                              sex=as.character(NA),
                                                                                                              ...){
                                                                                                           
Raise_Lgth(csObject,clObject,dbeOutput,spp=spp,taxon=taxon,sex=sex)

})














             