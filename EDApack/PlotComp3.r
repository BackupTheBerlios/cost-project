########################################################################################################################
##                                                                                                                    ##
## Plots of [volume of landings/effort] and [number of samples/ number of fish measured] by time and technical strata ##
##                                                                                                                    ##
##                          MM 14/01/2008                                                                             ##
########################################################################################################################

##save 'sole2' & 'GP' .RData files and change the path
setwd("C:/")

require (COSTcore)
load("sole2.RData") 



setGeneric("SampComp.plot", function(object1,object2,...){
	standardGeneric("SampComp.plot")
})

                                                                     
setMethod("SampComp.plot", signature(object1="csData",object2="clData"), function(object1,object2,TechStrat="commCat",TimeStrat="quarter",...){
#Tests on stratification parameters
if (!TechStrat%in%c("commCat","foCatNat","foCatEu5","foCatEu6")) stop("Wrong 'TechStrat' parameter!")
if (!TimeStrat%in%c("year","quarter","month")) stop("Wrong 'TimeStrat' parameter!")            #TimeStrat field is supposed to be numerical (for ordering)
#If TimeStrat="quarter" or "month", field must be put in HH
HHtime <- object1@hh$year
HHmonth <- as.numeric(sapply(object1@hh$date,function(x) strsplit(as.character(x),"-")[[1]][2]))
if (TimeStrat=="month") HHtime <- HHmonth
if (TimeStrat=="quarter") HHtime <- floor((HHmonth-0.1)/3)+1
object1@hh$HHtime <- HHtime
#insert 'Stratas' from HH to HL 
tab <- merge(object1@hl,object1@hh,all.x=TRUE)  
#total weight in the stratification
pdsPopTechTim <- tapply(object2@cl$landWt,list(object2@cl[,TechStrat],object2@cl[,TimeStrat]),sum,na.rm=TRUE)  
#NA modality and empty lines/columns taken off
pdsPopTechTim <- pdsPopTechTim[!dimnames(pdsPopTechTim)[[1]]%in%c("NA","<NA>"),!dimnames(pdsPopTechTim)[[2]]%in%c("NA","<NA>")]
pdsPopTechTim <- pdsPopTechTim[apply(pdsPopTechTim,1,function(x) !all(is.na(x))),apply(pdsPopTechTim,2,function(x) !all(is.na(x)))]  
pdsPopTechTim[is.na(pdsPopTechTim)] <- 0
#number of fish measured (for the same modalities as pdsPopTechTim's)
nbMesTechTim <- tapply(tab$lenNum,list(factor(tab[,TechStrat],levels=dimnames(pdsPopTechTim)[[1]]),factor(tab$HHtime,levels=dimnames(pdsPopTechTim)[[2]])),sum,na.rm=TRUE)
nbMesTechTim[is.na(nbMesTechTim)] <- 0
#Warnings on stratas modalities included in CS but not in CL (ex: commCat=51)
TS <- unique(as.character(tab[,TechStrat])) ; TS <- TS[!is.na(TS)] ; TS <- TS[!TS=="NA"] ; TS <- TS[!TS%in%dimnames(pdsPopTechTim)[[1]]]
TT <- unique(as.character(tab$HHtime)) ; TT <- TT[!is.na(TT)] ; TT <- TT[!TT=="NA"] ; TT <- TT[!TT%in%dimnames(pdsPopTechTim)[[2]]]
if (length(TS)!=0) warning(paste("Some CS technical strata modalities are not in CL table!!",paste(TS,collapse=",")))
if (length(TT)!=0) warning(paste("Some CS time strata modalities are not in CL table!!",paste(TT,collapse=",")))

TechMes <- apply(nbMesTechTim,1,sum,na.rm=TRUE)/sum(nbMesTechTim,na.rm=TRUE) ; TechPop <- apply(pdsPopTechTim,1,sum,na.rm=TRUE)/sum(pdsPopTechTim,na.rm=TRUE)
TimMes <- apply(nbMesTechTim,2,sum,na.rm=TRUE)/sum(nbMesTechTim,na.rm=TRUE) ; TimPop <- apply(pdsPopTechTim,2,sum,na.rm=TRUE)/sum(pdsPopTechTim,na.rm=TRUE)
dimVec <- c(length(TechMes),length(TimMes))
TechDf <- data.frame(str=factor(c(names(TechMes),names(TimMes)),levels=c(names(TimMes),names(TechMes))),mes=c(TechMes,TimMes),pop=c(TechPop,TimPop),grp=rep(1:2,dimVec),GRP=rep(c("Technical","Time"),dimVec))
mes <- TechDf$mes ; names(mes) <- TechDf$grp

#graphs 
require(lattice)
load("GP.RData")                                                                                                                               #<<<<---- to be replaced by 'data(...)'
dots <- list(...) ; if (is.null(dots$p.col)) dots$p.col <- "black" ; if (is.null(dots$p.bg)) dots$p.bg <- "lightblue" ; if (is.null(dots$l.col)) dots$l.col <- "red"
if (is.null(dots$lwd)) dots$lwd <- 2 ; if (is.null(dots$pch)) dots$pch <- 20    #graphical default values

sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "" ; if (is.null(dots$ylab)) dots$ylab <- "Frequency" 
if (is.null(dots$main)) dots$main <- "Relative Rates of Total Weight of Landings and Number of Fish Measured\nby Time and Technical Strata"

barchart(pop~str|GRP,data=TechDf,scales=list(x=list(relation="free"),font=dots$font.axis),main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),
         prepanel=function(x,y,...){
             x <- x[,drop=TRUE]
             prepanel.default.bwplot(x,y,...)},
         panel = function(x,y,mis=mes,...){
             x <- x[,drop=TRUE]
             panel.barchart(x,y,col=dots$p.bg,fill=dots$p.bg,...)
             panel.lines(type="o",mis[names(mis)==as.character(packet.number())],col=dots$l.col,lwd=dots$lwd,pch=dots$pch,cex=dots$cex,lty=dots$lty)},
         layout=1:2,ylim=c(0,1)) 
})





###########
##Example##
###########

#only market sampling data has 'commCat' values --> no need to subset for sampType=="M"
load("sole2.RData")                                                                                                                            #<<<<---- to be replaced by 'data(...)'
SampComp.plot(sole2.cs,sole2.cl)
SampComp.plot(sole2.cs,sole2.cl,TimeStrat="month")
SampComp.plot(sole2.cs,sole2.cl,l.col="steelblue",p.bg="gold",pch=15,lwd=1,lty=2)




object1 <- sole2.cs
object2 <- sole2.cl
#subset of object1 to quarters 2 & 4, and no 30 CC 
#on extrait les mois de hh
hhMonth <- sapply(object1@hh$date,function(x) as.numeric(strsplit(as.character(x),"-")[[1]][2])) 
object1@hh <- object1@hh[hhMonth%in%c(4:6,10:12),]
object1@sl <- object1@sl[object1@sl$trpCode%in%object1@hh$trpCode,] ; object1@hl <- object1@hl[object1@hl$trpCode%in%object1@hh$trpCode,]
object1@sl <- object1@sl[!object1@sl$commCat==30,] ; object1@hl <- object1@hl[!object1@hl$commCat==30,]

SampComp.plot(object1,object2)

