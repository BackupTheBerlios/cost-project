########################################################################################################################
##                                                                                                                    ##
## Plots of [volume of landings/effort] and [number of samples/ number of fish measured] by time and technical strata ##
##                                                                                                                    ##
##                          MM 29/01/2008                                                                             ##
########################################################################################################################

##save 'Sole' & 'GraphsPar' .RData files and change the path
setwd("C:/")

require (COSTcore)


setGeneric("SampComp.plot", function(object1,object2,...){
	standardGeneric("SampComp.plot")
})

                                                                     
setMethod("SampComp.plot", signature(object1="csData",object2="clData"), function(object1,object2,TechStrat="commCat",TimeStrat="quarter",SpaceStrat="area",show="all",...){
if ("all"%in%show) show <- c("Space","Technical","Time")
#Tests on stratification parameters
if (!TechStrat%in%c("commCat","foCatNat","foCatEu5","foCatEu6")) stop("Wrong 'TechStrat' parameter!")
if (!TimeStrat%in%c("year","quarter","month")) stop("Wrong 'TimeStrat' parameter!")            #TimeStrat field is supposed to be numerical (for ordering)
if (!SpaceStrat%in%c("area","rect")) stop("Wrong 'SpaceStrat' parameter!")           

#If TimeStrat="quarter" or "month", field must be put in HH
HHtime <- object1@hh$year
HHmonth <- as.numeric(sapply(object1@hh$date,function(x) strsplit(as.character(x),"-")[[1]][2]))
if (TimeStrat=="month") HHtime <- HHmonth
if (TimeStrat=="quarter") HHtime <- floor((HHmonth-0.1)/3)+1
object1@hh$HHtime <- HHtime
#insert 'Stratas' from HH to HL 
tab <- merge(object1@hl,object1@hh,all.x=TRUE)  
#total weight in the stratification
pdsPopTechTim <- tapply(object2@cl$landWt,list(object2@cl[,TechStrat],object2@cl[,TimeStrat],object2@cl[,SpaceStrat]),sum,na.rm=TRUE)  
#NA modality and empty lines/columns taken off
pdsPopTechTim <- pdsPopTechTim[!dimnames(pdsPopTechTim)[[1]]%in%c("NA","<NA>"),!dimnames(pdsPopTechTim)[[2]]%in%c("NA","<NA>"),!dimnames(pdsPopTechTim)[[3]]%in%c("NA","<NA>")]
pdsPopTechTim <- pdsPopTechTim[apply(pdsPopTechTim,1,function(x) !all(is.na(x))),apply(pdsPopTechTim,2,function(x) !all(is.na(x))),apply(pdsPopTechTim,3,function(x) !all(is.na(x)))]  
pdsPopTechTim[is.na(pdsPopTechTim)] <- 0
#number of fish measured (for the same modalities as pdsPopTechTim's)
nbMesTechTim <- tapply(tab$lenNum,list(factor(tab[,TechStrat],levels=dimnames(pdsPopTechTim)[[1]]),factor(tab$HHtime,levels=dimnames(pdsPopTechTim)[[2]]),
                                       factor(tab[,SpaceStrat],levels=dimnames(pdsPopTechTim)[[3]])),sum,na.rm=TRUE)
nbMesTechTim[is.na(nbMesTechTim)] <- 0

TechMes <- apply(nbMesTechTim,1,sum,na.rm=TRUE)/sum(nbMesTechTim,na.rm=TRUE) ; TechPop <- apply(pdsPopTechTim,1,sum,na.rm=TRUE)/sum(pdsPopTechTim,na.rm=TRUE)
TimMes <- apply(nbMesTechTim,2,sum,na.rm=TRUE)/sum(nbMesTechTim,na.rm=TRUE) ; TimPop <- apply(pdsPopTechTim,2,sum,na.rm=TRUE)/sum(pdsPopTechTim,na.rm=TRUE)
SpacMes <- apply(nbMesTechTim,3,sum,na.rm=TRUE)/sum(nbMesTechTim,na.rm=TRUE) ; SpacPop <- apply(pdsPopTechTim,3,sum,na.rm=TRUE)/sum(pdsPopTechTim,na.rm=TRUE)

dimVec <- c(length(SpacMes),length(TechMes),length(TimMes))
TechDf <- data.frame(str=factor(c(names(SpacMes),names(TechMes),names(TimMes)),levels=c(names(SpacMes),names(TimMes),names(TechMes))),
                     mes=c(SpacMes,TechMes,TimMes),pop=c(SpacPop,TechPop,TimPop),GRP=rep(c("Space","Technical","Time"),dimVec))
TechDf <- TechDf[TechDf$GRP%in%show,] ; TechDf$grp <- factor(TechDf$GRP,labels=1:length(unique(TechDf$GRP)))
mes <- TechDf$mes ; names(mes) <- TechDf$grp

#graphs 
require(lattice)
load("GraphsPar.RData")                                                                                                                               #<<<<---- to be replaced by 'data(...)'
dots <- list(...) ; if (is.null(dots$p.col)) dots$p.col <- "black" ; if (is.null(dots$p.bg)) dots$p.bg <- "lightblue" ; if (is.null(dots$l.col)) dots$l.col <- "red"
if (is.null(dots$lwd)) dots$lwd <- 2 ; if (is.null(dots$pch)) dots$pch <- 20    #graphical default values

sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "" ; if (is.null(dots$ylab)) dots$ylab <- "Frequency" 
if (is.null(dots$main)) dots$main <- paste("Relative Rates of Total Weight of Landings and Number of Fish Measured\nby",paste(unique(TechDf$GRP),collapse=", "),"Strata")

barchart(pop~str|GRP,data=TechDf,scales=list(x=list(relation="free"),font=dots$font.axis),main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),
         prepanel=function(x,y,...){
             x <- x[,drop=TRUE]
             prepanel.default.bwplot(x,y,...)},
         panel = function(x,y,mis=mes,...){
             x <- x[,drop=TRUE]
             panel.barchart(x,y,col=dots$p.bg,fill=dots$p.bg,...)
             panel.lines(type="o",mis[names(mis)==as.character(packet.number())],col=dots$l.col,lwd=dots$lwd,pch=dots$pch,cex=dots$cex,lty=dots$lty)},
         layout=c(1,length(unique(TechDf$GRP))),ylim=c(0,1)) 
})





###########
##Example##
###########

#only market sampling data has 'commCat' values --> no need to subset for sampType=="M"
load("Sole.RData")                                                                                                                            #<<<<---- to be replaced by 'data(...)'
SampComp.plot(sole3.cs,sole3.cl)
SampComp.plot(sole3.cs,sole3.cl,TimeStrat="month",show=c("Time","Technical"))
SampComp.plot(sole3.cs,sole3.cl,l.col="steelblue",p.bg="gold",pch=15,lwd=1,lty=2,show=c("Time","Technical"))



