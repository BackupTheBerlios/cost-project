########################################################
##                                                    ##
## Plots of Volume of Landings/Discards per Haul/Trip ##
##                                                    ##
##       MM 14/01/2008                                ##
########################################################

##save 'sole2' & 'GP' .RData files and change the path

setwd("C:/")
require(COSTcore)


setClass("LD.Vol",representation(fraction="character",species="character",TechStrat="character",TimeStrat="character",VolFO_FDTR="list",MeanFO_FDTR="numeric",VolFD_TR="list",
                  MeanFD_TR="numeric",StratDF="data.frame"),
                  prototype(fraction=character(),species=character(),TechStrat=character(),TimeStrat=character(),VolFO_FDTR=list(),MeanFO_FDTR=numeric(),
                  VolFD_TR=list(),MeanFD_TR=numeric(),StratDF=data.frame()))

                     


setGeneric("LD.Volume", function(object,...){
	standardGeneric("LD.Volume")
})




setMethod("LD.Volume", signature(object="csData"), function(object,fraction="all",species="all",TechStrat="gear",TimeStrat="quarter",...){   #TimeStrat = "year", "quarter", "month"   
                                                                                                                            #TechStrat = "gear", "foCatNat", "foCatEu5", "foCatEu6"
op.sub <- object@hh[object@hh$sampType=="S",]                                     
capt.sub <- object@sl[paste(object@sl$trpCode,object@sl$staNum,sep="::")%in%paste(op.sub$trpCode,op.sub$staNum,sep="::"),]  

if (!"all"%in%fraction) capt.sub <- capt.sub[capt.sub$catchCat%in%fraction,]                                                                                                                                              
if (!"all"%in%species) capt.sub <- capt.sub[capt.sub$spp%in%species,]                                                                                                                                              
                                                                                                                                                                  
op.sub$trpCode <- factor(op.sub$trpCode) ; op.sub$date <- factor(op.sub$date) ; op.sub$staNum <- factor(op.sub$staNum)

#If TimeStrat="quarter" or "month", field must be put in HH
HHtime <- op.sub$year
HHmonth <- as.numeric(sapply(op.sub$date,function(x) strsplit(as.character(x),"-")[[1]][2]))
if (TimeStrat=="month") HHtime <- HHmonth
if (TimeStrat=="quarter") HHtime <- floor((HHmonth-0.1)/3)+1
op.sub$HHtime <- HHtime

#on determine une valeur technique et temporelle par maree (l'occurence qu'on retrouve le plus souvent)
TechVal <- tapply(op.sub[,TechStrat],list(op.sub$trpCode),function(x) names(sort(table(x),decreasing=TRUE)[1]))
TimeVal <- tapply(op.sub$HHtime,list(op.sub$trpCode),function(x) names(sort(table(x),decreasing=TRUE)[1]))

#Number of sampled fishing days by trip
tablEch <- op.sub[op.sub$foVal=="F",]                                                  
d_i <- tapply(tablEch$date,list(tablEch$trpCode),function(x) length(unique(x)))                 
SampTest <- !is.na(d_i)                                                 

#Number of FOs by fishing day and by trip 
M_ik <- tapply(op.sub$staNum,list(op.sub$trpCode,op.sub$date),function(x) length(unique(x)))

#Number of sampled FOs by fishing day and by trip
m_ik <- tapply(tablEch$staNum,list(tablEch$trpCode,tablEch$date),function(x) length(unique(x)))          


tabl1 <- merge(tablEch,aggregate(capt.sub$wt,list(trpCode=capt.sub$trpCode,staNum=capt.sub$staNum),sum),all.x=TRUE) #essentially to homogenize vectors sizes                
names(tabl1)[ncol(tabl1)] <- "wt"    
tabl1$wt[is.na(tabl1$wt)] <- 0             #species not present in fraction of the FO                                    
  
y_ikj <- tapply(tabl1$wt,list(tabl1$trpCode,tabl1$date,tabl1$staNum),sum,na.rm=TRUE)
  
                                                                                        
y_ikj_hat <- split(tabl1$wt,paste(tabl1$trpCode,tabl1$date,sep="::"),drop=TRUE)
y_ik <- unlist(lapply(y_ikj_hat,mean))
y_IK <- M_ik*apply(y_ikj,c(1,2),sum,na.rm=TRUE)/m_ik
mat <- matrix(rep(dimnames(y_IK)[[1]],ncol(y_IK)),ncol=ncol(y_IK)) ; mat[is.na(y_IK)] <- NA
y_ik_hat <- split(y_IK,mat,drop=TRUE) 
y_i <- unlist(lapply(y_ik_hat,mean))

StratDF <- data.frame(trp=names(y_i),tech=TechVal[names(y_i)],tim=TimeVal[names(y_i)])

invisible(new("LD.Vol",fraction=fraction,species=species,TechStrat=TechStrat,TimeStrat=TimeStrat,VolFO_FDTR=y_ikj_hat,MeanFO_FDTR=y_ik,VolFD_TR=y_ik_hat,MeanFD_TR=y_i,StratDF=StratDF))
})  




################################################################################
################################################################################


setGeneric("plot.FD", function(x,...){
	standardGeneric("plot.FD")
})


setMethod("plot.FD", signature(x="LD.Vol"), function(x,strat=FALSE,...){       
require(lattice) 
load("GP.RData")                                                                                                           #<<<<---- to be replaced by 'data(...)'
dots <- list(...)
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Trip Code" ; if (is.null(dots$ylab)) dots$ylab <- "Weight (Kg)" 
if ((is.null(dots$main))&(!strat)) dots$main <- paste("Mean Weight by Fishing Day for each Trip\nSpecies :",paste(x@species,collapse=", "),"    Fraction :",paste(x@fraction,collapse=", "))
if ((is.null(dots$main))&(strat)) dots$main <- paste("Mean Weight by Fishing Day for each Trip, by Time and Technical strata\nSpecies :",
                                                      paste(x@species,collapse=", "),"    Fraction :",paste(x@fraction,collapse=", "),"\nTime Strata :",
                                                      x@TimeStrat,"    Technical Strata :",x@TechStrat)
 
df <- cbind.data.frame(x@StratDF,data.frame(bb=x@MeanFD_TR))
TECH <- levels(factor(df$tech))
if (!strat) {
xyplot(bb~trp,data=df,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),
                     scales=list(font=dots$font.axis,x=list(rot=dots$rot)),
                     pch=dots$pch[1],fill=dots$bg,cex=dots$p.cex[1],col=dots$col[1])
} else {                     
xyplot(bb~trp|tim,data=df,groups=tech,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),
                     scales=list(font=dots$font.axis,x=list(relation="free",rot=dots$rot)),key=list(points=list(pch=dots$pch[1],fill=dots$p.bg[1:length(TECH)],cex=dots$p.cex[1],col=dots$col[1]),
              text=list(TECH,font=dots$font.lab),space="right",columns=1,border=TRUE),
              prepanel=function(x,y,...){
                  x <- x[,drop=TRUE]
                  prepanel.default.xyplot(x,y,...)},
              panel = function(x,y,...){
                  x <- x[,drop=TRUE]
                  panel.xyplot(x,y,pch=dots$pch[1],fill=dots$p.bg,cex=dots$p.cex[1],col=dots$col[1],...)})
}
})



################################################################################
################################################################################


setGeneric("boxplot.FD", function(x,...){
	standardGeneric("boxplot.FD")
})


setMethod("boxplot.FD", signature(x="LD.Vol"), function(x,...){       
  
load("C:/GP.RData")                                                                                                             #<<<<---- to be replaced by 'data(...)'
dots <- list(...) ; if (is.null(dots$pch)) dots$pch <- 19
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Trip Code" ; if (is.null(dots$ylab)) dots$ylab <- "Weight (Kg)" ; if (is.null(dots$main)) dots$main <- "Weight by Fishing Day for each Trip"

vec <- unlist(x@VolFD_TR) ; nvec <- unlist(lapply(x@VolFD_TR,length)) 
df <- data.frame(aa=rep(names(nvec),nvec),bb=vec)

bwplot(bb~aa,data=df,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),
                     pch=dots$pch[1],fill=dots$p.bg[1],scales=list(font=dots$font.axis,x=list(rot=dots$rot)),
                     par.settings=list(box.rectangle=list(col=dots$col[1]),box.umbrella=list(col=dots$col[1],lty=dots$lty[1]),
                     plot.symbol=list(col=dots$col[1])))

})



################################################################################
################################################################################


setGeneric("plot.FO", function(x,...){
	standardGeneric("plot.FO")
})


setMethod("plot.FO", signature(x="LD.Vol"), function(x,...){       
require(lattice)  
load("GP.RData")                                                                                                                 ##<<<<---- to be replaced by 'data(...)'
dots <- list(...) ; if (is.null(dots$rot)) dots$rot <- 0
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Fishing day" ; if (is.null(dots$ylab)) dots$ylab <- "Weight (Kg)" 
if (is.null(dots$main)) dots$main <- paste("Mean Weight by FO for each Fishing day and Trip\nSpecies :",paste(x@species,collapse=", "),"    Fraction :",paste(x@fraction,collapse=", "))
 
mat <- t(sapply(names(x@MeanFO_FDTR),function(x) strsplit(x,"::")[[1]]))

df <- data.frame(trpCode=mat[,1],date=mat[,2],Fday=as.numeric(unlist(tapply(mat[,1],list(mat[,1]),function(x) 1:length(x)))),val=as.numeric(x@MeanFO_FDTR))
   
xyplot(val~Fday|trpCode,data=df,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),
                     scales=list(font=dots$font.axis,x=list(rot=dots$rot)),par.strip.text=list(font=dots$font.lab),
                     pch=dots$pch[1],fill=dots$bg,cex=dots$p.cex[1],col=dots$col[1])
})




################################################################################
################################################################################


setGeneric("boxplot.FO", function(x,...){
	standardGeneric("boxplot.FO")
})


setMethod("boxplot.FO", signature(x="LD.Vol"), function(x,...){       
require(lattice)  
load("GP.RData")                                                                                                                 ##<<<<---- to be replaced by 'data(...)'
dots <- list(...) ; if (is.null(dots$rot)) dots$rot <- 0 ; if (is.null(dots$pch)) dots$pch <- 19
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Fishing day" ; if (is.null(dots$ylab)) dots$ylab <- "Weight (Kg)" 
if (is.null(dots$main)) dots$main <- paste("Weight by Fishing Operation for each Fishing day and Trip\nSpecies :",paste(x@species,collapse=", "),"    Fraction :",paste(x@fraction,collapse=", "))

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






###########
# Example # 
###########



load("sole2.RData")
object <- sole2.cs
#only sea sampling data is kept
object@tr <- object@tr[object@tr$sampType=="S",] 
object@hh <- object@hh[object@hh$sampType=="S",] 
object@sl <- object@sl[object@sl$sampType=="S",]
object@hl <- object@hl[object@hl$sampType=="S",]


x <- LD.Volume(object,fraction="LAN",species="SOL")

plot.FD(x,col="blue",bg="steelblue",p.cex=1.4)
plot.FD(x,strat=TRUE)

boxplot.FD(x)
   
plot.FO(x,col="blue",bg="steelblue",p.cex=1.4)

boxplot.FO(x)
 

 