##################################################################
##                                                              ##
## Plots of Volume of Landings/Discards per FO/Fishing Day/Trip ##
##                                                              ##
##                      MM 07/02/2008                           ##
##################################################################


setClass("LD.Vol",representation(fraction="character",species="character",Strata="list",VolFO_FDTR="list",MeanFO_FDTR="numeric",
                  VolFD_TR="list",MeanFD_TR="numeric"),
                  prototype(fraction=character(),species=character(),Strata=list(),VolFO_FDTR=list(),MeanFO_FDTR=numeric(),
                  VolFD_TR=list(),MeanFD_TR=numeric()))


################################################################################
################################################################################
                     

setGeneric("LD.Volume", function(object,                                   
                                 species,
                                 fraction="LAN",
                                 TimeStrat=NULL,
                                 TechStrat=NULL,
                                 SpaceStrat=NULL,...){
	standardGeneric("LD.Volume")
})

    
                                                       


setMethod("LD.Volume", signature(object="csData"), function(object,
                                                            species,         #une espèce seulement ou "all" si la table SL n'en contient qu'une
                                                            fraction="LAN",
                                                            TimeStrat=NULL,
                                                            TechStrat=NULL,
                                                            SpaceStrat=NULL,...){  


                                                          
op.sub <- object@hh[object@hh$sampType=="S",]   #only sea sampling data    
                                                                    
capt.sub <- object@sl[paste(object@sl$trpCode,object@sl$staNum,sep="::")%in%paste(op.sub$trpCode,op.sub$staNum,sep="::"),]  

un <- unique(as.character(object@sl$spp)) ; un <- un[!is.na(un)]                                                                    
if (species=="all") {if (length(un)>1) {warning("Several species in SL table!! Only the first one will be taken into account!")}
                     species <- un[1]} 

capt.sub <- capt.sub[capt.sub$catchCat%in%fraction,]                                                                                                                                              
capt.sub <- capt.sub[capt.sub$spp%in%species,]                                                                                                                                              
                                                                                                                                                                  
op.sub$trpCode <- factor(op.sub$trpCode) ; op.sub$date <- factor(op.sub$date) ; op.sub$staNum <- factor(op.sub$staNum)

#If TimeStrat="semester", "quarter" or "month", field must be put in HH
if (!is.null(TimeStrat)) {
HHmonth <- as.numeric(sapply(op.sub$date,function(x) strsplit(as.character(x),"-")[[1]][2]))
if (TimeStrat=="month") op.sub$month <- HHmonth
if (TimeStrat=="quarter") op.sub$quarter <- ceiling(HHmonth/3)
if (TimeStrat=="semester") op.sub$semester <- ceiling(HHmonth/6) 
}

#les champs de stratification dans op.sub sont transformés en facteurs
if (!is.null(TimeStrat)) op.sub[,TimeStrat] <- factor(op.sub[,TimeStrat])
if (!is.null(TechStrat)) op.sub[,TechStrat] <- factor(op.sub[,TechStrat])
if (!is.null(SpaceStrat)) op.sub[,SpaceStrat] <- factor(op.sub[,SpaceStrat])
 


#Number of sampled fishing days by trip, tech,time,space
#tablEch <- op.sub[op.sub$foVal=="V",]                                                  
expr1 <- paste(",tabSamp$",c(TimeStrat,TechStrat,SpaceStrat),sep="",collapse="") ; if (expr1==",tabSamp$") expr1 <- ""                 
expr2 <- paste(",op.sub$",c(TimeStrat,TechStrat,SpaceStrat),sep="",collapse="") ; if (expr2==",op.sub$") expr2 <- ""                 
expr3 <- paste(",tabl1$",c(TimeStrat,TechStrat,SpaceStrat),sep="",collapse="") ; if (expr3==",tabl1$") expr3 <- ""                 


#eval(parse('',text=paste("d_i <- tapply(tablEch$date,list(tablEch$trpCode",expr1,"),function(x) length(unique(x)))",sep="")))                 
#SampTest <- !is.na(d_i)                                                 
#
#Number of FOs by fishing day, by trip, by tech,time,space 
eval(parse('',text=paste("M_ik <- tapply(op.sub$staNum,list(op.sub$trpCode,op.sub$date",expr2,"),function(x) length(unique(x)))",sep="")))

#Number of sampled FOs by fishing day, by trip, by tech,time,space 
#Fo is considered sampled if foVal=="V" and : catchReg=="All"  &&  sppReg=="All"
#                                             catchReg==fraction && sppReg=="All"
#                                             catchReg=="All"  &&  sppReg=="Par"  && species in SL table for the FO and the fraction
#                                             catchReg==fraction  &&  sppReg=="Par"  && species in SL table for the FO and the fraction 
if (fraction=="LAN") fract <- "Lan" else fract <- "Dis" 
capt.sub$ind <- 1  
tabSamp <- merge(op.sub,capt.sub[,c("sampType","landCtry","vslFlgCtry","year","proj","trpCode","staNum","ind")],all.x=TRUE)
tabSamp$ind[is.na(tabSamp$ind)] <- 0 ; tabSamp$ind2 <- tabSamp$ind   #sampling indicator
tabSamp$ind2[tabSamp$catReg=="All" & tabSamp$sppReg=="All"] <- 1
tabSamp$ind2[tabSamp$catReg==fract & tabSamp$sppReg=="All"] <- 1
tabSamp$ind2[tabSamp$catReg=="All" & tabSamp$sppReg=="Par" & tabSamp$ind==1] <- 1
tabSamp$ind2[tabSamp$catReg==fract & tabSamp$sppReg=="Par" & tabSamp$ind==1] <- 1
tabSamp$ind2[tabSamp$foVal!="V"] <- 0

eval(parse('',text=paste("m_ik <- tapply(tabSamp$ind2,list(tabSamp$trpCode,tabSamp$date",expr1,"),sum)",sep="")))          


tabl1 <- merge(tabSamp,aggregate(capt.sub$wt,list(sampType=capt.sub$sampType,landCtry=capt.sub$landCtry,vslFlgCtry=capt.sub$vslFlgCtry,year=capt.sub$year,
                                                  proj=capt.sub$proj,trpCode=capt.sub$trpCode,staNum=capt.sub$staNum),sum),all.x=TRUE) #essentially to homogenize vectors sizes                
names(tabl1)[ncol(tabl1)] <- "wt"    
tabl1$wt[is.na(tabl1$wt)] <- 0                                                
  
#on se refère à la colonne ind2 pour ne garder de tabl1 que ce qui est considéré comme échantillonné  
tabl1 <- tabl1[tabl1$ind2==1,]  
  
eval(parse('',text=paste("y_ikj <- tapply(tabl1$wt,list(tabl1$trpCode,tabl1$date",expr3,"),sum,na.rm=TRUE)",sep="")))
  
                                                                                        
y_ikj_hat <- split(tabl1$wt,paste(tabl1$trpCode,tabl1$date,sep="::"),drop=TRUE)
y_ik <- unlist(lapply(y_ikj_hat,mean))

y_IK <- M_ik*y_ikj/m_ik
ll <- sum(!is.null(TechStrat),!is.null(TimeStrat),!is.null(SpaceStrat))
indic <- ll+2 ; val <- expand.grid(dimnames(y_IK)[1:indic]) ; valChar <- apply(val[,-2,drop=FALSE],1,function(x) paste(as.character(x),collapse="::"))
MAT <- array(valChar,dim=dim(y_IK)) ; MAT[is.na(y_IK)] <- NA

y_ik_hat <- split(y_IK,MAT,drop=TRUE) 
y_i <- unlist(lapply(y_ik_hat,mean))


invisible(new("LD.Vol",fraction=fraction,species=species,Strata=list(TimeStrat=TimeStrat,TechStrat=TechStrat,SpaceStrat=SpaceStrat),
               VolFO_FDTR=y_ikj_hat,MeanFO_FDTR=y_ik,VolFD_TR=y_ik_hat,MeanFD_TR=y_i))
})  


################################################################################
################################################################################


setGeneric("FD.plot", function(x,
                               groups=NULL,...){
	standardGeneric("FD.plot")
})


setMethod("FD.plot", signature(x="LD.Vol"), function(x,
                                                     groups=NULL,...){       

stratas <- c("TimeStrat","TechStrat","SpaceStrat")
TimeStrat <- x@Strata$TimeStrat ; TechStrat <- x@Strata$TechStrat ; SpaceStrat <- x@Strata$SpaceStrat 
index <- c(TimeStrat,TechStrat,SpaceStrat)
data(GraphsPar)                                                                                                       
dots <- list(...)
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Trip Code" ; if (is.null(dots$ylab)) dots$ylab <- "Weight (Kg)" 


if (all(is.null(TimeStrat),is.null(TechStrat),is.null(SpaceStrat))) {

if (is.null(dots$main)) dots$main <- paste("Mean Weight by Fishing Day for each Trip\nSpecies :",paste(x@species,collapse=", "),"    Fraction :",paste(x@fraction,collapse=", "))
df <- data.frame(trp=names(x@MeanFD_TR),bb=as.numeric(x@MeanFD_TR))

xyplot(bb~trp,data=df,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),
       scales=list(font=dots$font.axis,x=list(rot=dots$rot)),pch=dots$pch[1],fill=dots$p.bg[1],cex=dots$p.cex[1],col=dots$col[1])
                     
} else {    

if (is.null(dots$main)) dots$main <- paste("Mean Weight by Fishing Day for each Trip, by Time and Technical strata\nSpecies :",
                                           paste(x@species,collapse=", "),"    Fraction :",paste(x@fraction,collapse=", "),"\n",
                                           paste("Time Strata :",x@Strata$TimeStrat)[!is.null(x@Strata$TimeStrat)],
                                           paste("   Technical Strata :",x@Strata$TechStrat)[!is.null(x@Strata$TechStrat)],
                                           paste("   Space Strata :",x@Strata$SpaceStrat)[!is.null(x@Strata$SpaceStrat)])
datas <-  x@MeanFD_TR
df <- as.data.frame(do.call("rbind",lapply(names(datas),function(x) strsplit(x,"::")[[1]])))
names(df) <- c("trp",x@Strata$TimeStrat,x@Strata$TechStrat,x@Strata$SpaceStrat)
df$bb <- as.numeric(x@MeanFD_TR)

strip.col <- trellis.par.get("strip.background")$col
      
if (is.null(groups)) {
          
eval(parse('',text=paste("xyplot(bb~trp|",paste(index,collapse="*"),",data=df,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),",
                         "ylab=list(dots$ylab,font=dots$font.lab),par.strip.text=list(font=dots$font.lab),scales=list(font=dots$font.axis,x=list(relation=\"free\",rot=dots$rot)),",
                         "key=list(points=list(pch=15,cex=dots$p.cex[1],col=strip.col[1:length(index)]),text=list(index),font=dots$font.lab,columns=1,border=TRUE,space=\"right\"),",
                         "prepanel=function(x,y,...){x <- x[,drop=TRUE] ; prepanel.default.xyplot(x,y,...)},",
                         "panel = function(x,y,...){x <- x[,drop=TRUE] ; panel.xyplot(x,y,pch=dots$pch[1],fill=dots$p.bg[1],cex=dots$p.cex[1],col=dots$col[1],...)})",sep="")))

} else {

indexStr <- index[!index%in%eval(parse('',text=groups))] ; l1 <- length(indexStr)
LEV <- levels(df[,eval(parse('',text=groups))]) ; l2 <- length(LEV)
groups <- eval(parse('',text=groups))

eval(parse('',text=paste("xyplot(bb~trp",paste("|",paste(indexStr,collapse="*"),sep="")[l1>0],",data=df,groups=",groups,",main=list(dots$main,font=dots$font.main),",
                         "xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),scales=list(font=dots$font.axis,x=list(relation=\"free\",rot=dots$rot)),",
                         "key=list(points=list(pch=c(rep(dots$pch[1],l2),NA,",c("rep(15,l1)","NA")[c((l1>0),(l1==0))],"),fill=dots$p.bg[1:l2],cex=dots$p.cex[1],lwd=dots$p.lwd[1],",
                         "col=c(rep(dots$col[1],l2),NA",",strip.col[1:l1]"[l1>0],")),text=list(c(LEV,\"\",\"",paste(indexStr,collapse="\",\""),"\")),title=\"",groups,"\",",
                         "cex.title=0.8,space=\"right\",font=dots$font.lab,columns=1,border=TRUE),par.strip.text=list(font=dots$font.lab),",
                         "prepanel=function(x,y,...){x <- x[,drop=TRUE] ; prepanel.default.xyplot(x,y,...)},",
                         "panel = function(x,y,...){x <- x[,drop=TRUE] ; panel.xyplot(x,y,pch=dots$pch[1],fill=dots$p.bg","[1]"[l2==0],",cex=dots$p.cex[1],col=dots$col[1],...)})",sep="")))
}}
})


################################################################################
################################################################################


setGeneric("FD.boxplot", function(x,...){
	standardGeneric("FD.boxplot")
})


setMethod("FD.boxplot", signature(x="LD.Vol"), function(x,...){       
  
data(GraphsPar)                                                                                                        
dots <- list(...) ; if (is.null(dots$pch)) dots$pch <- 19
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Trip Code" ; if (is.null(dots$ylab)) dots$ylab <- "Weight (Kg)" 
if (is.null(dots$main)) dots$main <- paste("Weight by Fishing Day for each Trip\nSpecies :",
                                           paste(x@species,collapse=", "),"    Fraction :",paste(x@fraction,collapse=", "))

obj <- x@VolFD_TR ; names(obj) <- sapply(names(x@VolFD_TR),function(x) strsplit(x,"::")[[1]][1])
vec <- unlist(x@VolFD_TR) ; nvec <- unlist(lapply(x@VolFD_TR,length)) 
df <- data.frame(aa=rep(names(obj),nvec),bb=vec)

bwplot(bb~aa,data=df,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),
       pch=dots$pch[1],fill=dots$p.bg[1],scales=list(font=dots$font.axis,x=list(rot=dots$rot)),
       par.settings=list(box.rectangle=list(col=dots$col[1]),box.umbrella=list(col=dots$col[1],lty=dots$lty[1]),
       plot.symbol=list(col=dots$col[1])))

})


################################################################################
################################################################################


setGeneric("FO.plot", function(x,...){
	standardGeneric("FO.plot")
})


setMethod("FO.plot", signature(x="LD.Vol"), function(x,...){       

data(GraphsPar)                                                                                                   
dots <- list(...) ; if (is.null(dots$rot)) dots$rot <- 0
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Fishing day" ; if (is.null(dots$ylab)) dots$ylab <- "Weight (Kg)" 
if (is.null(dots$main)) dots$main <- paste("Mean Weight by FO for each Fishing day and Trip\nSpecies :",
                                           paste(x@species,collapse=", "),"    Fraction :",paste(x@fraction,collapse=", "))
 
mat <- t(sapply(names(x@MeanFO_FDTR),function(x) strsplit(x,"::")[[1]]))

df <- data.frame(trpCode=mat[,1],date=mat[,2],Fday=as.numeric(unlist(tapply(mat[,1],list(mat[,1]),function(x) 1:length(x)))),val=as.numeric(x@MeanFO_FDTR))
   
xyplot(val~Fday|trpCode,data=df,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),
       scales=list(font=dots$font.axis,x=list(rot=dots$rot)),par.strip.text=list(font=dots$font.lab),
       pch=dots$pch[1],fill=dots$p.bg[1],cex=dots$p.cex[1],col=dots$col[1])
})


################################################################################
################################################################################


setGeneric("FO.boxplot", function(x,...){
	standardGeneric("FO.boxplot")
})


setMethod("FO.boxplot", signature(x="LD.Vol"), function(x,...){       

data(GraphsPar)                                                                                                        
dots <- list(...) ; if (is.null(dots$rot)) dots$rot <- 0 ; if (is.null(dots$pch)) dots$pch <- 19
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Fishing day" ; if (is.null(dots$ylab)) dots$ylab <- "Weight (Kg)" 
if (is.null(dots$main)) dots$main <- paste("Weight by Fishing Operation for each Fishing day and Trip\nSpecies :",
                                           paste(x@species,collapse=", "),"    Fraction :",paste(x@fraction,collapse=", "))

vec <- unlist(x@VolFO_FDTR) ; nvec <- unlist(lapply(x@VolFO_FDTR,length)) ; nvec2 <- rep(names(nvec),nvec)
mat <- t(sapply(nvec2,function(x) strsplit(x,"::")[[1]]))

nvec3 <- sapply(unique(nvec2),function(x) strsplit(x,"::")[[1]][1])
ind1 <- as.numeric(unlist(tapply(nvec3,list(nvec3),function(x) 1:length(x))))
ind2 <- as.numeric(unlist(tapply(nvec2,list(nvec2),length)))

df <- data.frame(trpCode=mat[,1],date=mat[,2],Fday=rep(ind1,ind2),val=vec)
df$Fday <- as.factor(df$Fday)

bwplot(val~Fday|trpCode,data=df,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),
       pch=dots$pch[1],fill=dots$p.bg[1],scales=list(font=dots$font.axis,x=list(rot=dots$rot)),par.strip.text=list(font=dots$font.lab),
       par.settings=list(box.rectangle=list(col=dots$col[1]),box.umbrella=list(col=dots$col[1],lty=dots$lty[1]),
       plot.symbol=list(col=dots$col[1])))
})

