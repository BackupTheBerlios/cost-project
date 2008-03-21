#################################################################################################
##                                                                                             ##
## Plots of volume of landings and number of fish measured by time, technical and space strata ##
##                                                                                             ##
##                                   MM 06/02/2008                                             ##
#################################################################################################


setGeneric("SampComp.plot", function(object1,
                                     object2,
                                     var1="lenNum",
                                     var2="landWt",
                                     TimeStrat="quarter",
                                     TechStrat="commCat",
                                     SpaceStrat="area",
                                     show="all",
                                     ...){
	standardGeneric("SampComp.plot")
})

         
                                                                     
setMethod("SampComp.plot", signature(object1="csData",object2="clData"), function(object1,
                                                                                  object2,
                                                                                  var1="lenNum",
                                                                                  var2="landWt",
                                                                                  TimeStrat="quarter",
                                                                                  TechStrat="commCat",
                                                                                  SpaceStrat="area",
                                                                                  show="all",
                                                                                  ...){
if ("all"%in%show) show <- c("Space","Technical","Time")
TechInd <- is.null(TechStrat) ; TimeInd <- is.null(TimeStrat) ; SpaceInd <- is.null(SpaceStrat)
if (all(c(TechInd,TimeInd,SpaceInd))) stop("No stratification!!")
#Tests on stratification parameters
if (!TechInd) {if (!TechStrat%in%c("commCat","foCatNat","foCatEu5","foCatEu6")) stop("Wrong 'TechStrat' parameter!")}
if (!TimeInd) {if (!TimeStrat%in%c("year","quarter","month")) stop("Wrong 'TimeStrat' parameter!")}           
if (!SpaceInd) {if (!SpaceStrat%in%c("area","rect")) stop("Wrong 'SpaceStrat' parameter!")}           

#If TimeStrat="quarter" or "month", field must be put in HH
HHtime <- as.numeric(as.character(object1@hh$year))
HHmonth <- as.numeric(sapply(object1@hh$date,function(x) strsplit(as.character(x),"-")[[1]][2]))
if (!TimeInd) {if (TimeStrat=="month") HHtime <- HHmonth }
if (!TimeInd) {if (TimeStrat=="quarter") HHtime <- floor((HHmonth-0.1)/3)+1 }
HHtime <- factor(HHtime,exclude=NULL)  #sorting time occurrences
object1@hh$HHtime <- HHtime 
#insert 'Stratas' from HH to HL 
tab <- merge(object1@hl,object1@hh,all.x=TRUE)  


tabPop <- object2@cl


if (!TimeInd) tab[,"HHtime"] <- factor(tab[,"HHtime"],exclude="NA") 
if (!SpaceInd) tab[,SpaceStrat] <- factor(tab[,SpaceStrat],exclude="NA")
if (!TechInd) tab[,TechStrat] <- factor(tab[,TechStrat],exclude="NA") 
#population value in the stratification
eval(parse('',text=paste("pdsPopTechTim <- tapply(tabPop$",var2,",list(",
                         paste(c("factor(tabPop[,TimeStrat],exclude=NULL)"[!TimeInd],
                                 "factor(tabPop[,TechStrat],exclude=NULL)"[!TechInd],
                                 "factor(tabPop[,SpaceStrat],exclude=NULL)"[!SpaceInd]),collapse=","),"),sum,na.rm=TRUE)",sep=""))) 

index <- c(1:sum(!TimeInd,!TechInd,!SpaceInd))

eval(parse(text=paste("pdsPopTechTim <- pdsPopTechTim[",paste("apply(pdsPopTechTim,",
                          index,",function(x) !all(is.na(x)))",sep="",collapse=","),",drop=FALSE]",sep="")))  
                          
pdsPopTechTim[is.na(pdsPopTechTim)] <- 0

#sampled value in the stratification
mod <- c("\"HHtime\"","TechStrat","SpaceStrat")[c(!TimeInd,!TechInd,!SpaceInd)]
mod2 <- c("\"Time\"","\"Technical\"","\"Space\"")[c(!TimeInd,!TechInd,!SpaceInd)]

eval(parse('',text=paste("nbMesTechTim <- tapply(tab$lenNum,list(",paste("factor(tab[,",mod,"],levels=dimnames(pdsPopTechTim)[[",
                          index,"]],exclude=NULL)",sep="",collapse=","),"),sum,na.rm=TRUE)",sep="")))

nbMesTechTim[is.na(nbMesTechTim)] <- 0

eval(parse('',text=paste("Mes",index," <- apply(nbMesTechTim,",index,",sum,na.rm=TRUE)/sum(nbMesTechTim,na.rm=TRUE)",sep="",collapse=";")))
eval(parse('',text=paste("Pop",index," <- apply(pdsPopTechTim,",index,",sum,na.rm=TRUE)/sum(pdsPopTechTim,na.rm=TRUE)",sep="",collapse=";")))


eval(parse('',text=paste("Mes",index," <- Mes",index,"[!names(Mes",index,")%in%c(\"NA\",\"<NA>\",NA)]",sep="",collapse=";")))
eval(parse('',text=paste("Mes",index," <- Mes",index,"/sum(Mes",index,",na.rm=TRUE)",sep="",collapse=";")))
eval(parse('',text=paste("Pop",index," <- Pop",index,"[!names(Pop",index,")%in%c(\"NA\",\"<NA>\",NA)]",sep="",collapse=";")))
eval(parse('',text=paste("Pop",index," <- Pop",index,"/sum(Pop",index,",na.rm=TRUE)",sep="",collapse=";")))


eval(parse('',text=paste("dimVec <- c(",paste("length(Mes",index,")",sep="",collapse=","),")",sep="")))

int1 <- paste("names(Mes",index,")",sep="",collapse=",") ; int2 <- paste("Mes",index,sep="",collapse=",")
int3 <- paste("Pop",index,sep="",collapse=",") ; int4 <- paste(mod2,sep="",collapse=",") 

eval(parse('',text=paste("TechDf <- data.frame(str=factor(c(",int1,"),levels=c(",int1,")),mes=c(",int2,"),pop=c(",int3,"),GRP=rep(c(",int4,"),dimVec))",sep="")))

TechDf <- TechDf[TechDf$GRP%in%show,] ; TechDf$grp <- factor(TechDf$GRP,labels=1:length(unique(TechDf$GRP)))
TechDf$GRPplus <- factor(TechDf$GRP,levels=c("Space","Technical","Time"),labels=c(paste("Space =",SpaceStrat),paste("Technical =",TechStrat),paste("Time =",TimeStrat)))
mes <- TechDf$mes ; names(mes) <- TechDf$grp

#graph display 
data(GraphsPar)                                                                                                                           
dots <- list(...) 
if (is.null(dots$p.col)) dots$p.col <- "black" ; if (is.null(dots$p.bg)) dots$p.bg <- "lightblue" ; if (is.null(dots$l.col)) dots$l.col <- "red"
if (is.null(dots$lwd)) dots$lwd <- 2 ; if (is.null(dots$pch)) dots$pch <- 20  

sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "" ; if (is.null(dots$ylab)) dots$ylab <- "Frequency" 
if (is.null(dots$main)) dots$main <- paste("Relative Rates of \"",var1,"\" and \"",var2,"\" variables\nby ",paste(unique(TechDf$GRP),collapse=", "),"Strata",sep="")

barchart(pop~str|GRPplus,data=TechDf,scales=list(x=list(relation="free"),font=dots$font.axis),main=list(dots$main,font=dots$font.main),
         xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),par.strip.text=list(font=dots$font.lab),
         key =list(lines=list(pch=c(15,1),type=c("p","l"),col=c(dots$p.bg[1],dots$l.col[1]),lwd=c(2,dots$lwd[1]),cex=c(1.2,dots$cex[1]),lty=c(1,dots$lty[1])),
                   text=list(c(var2,var1)),font=dots$font.lab,space="right",columns=1,border=TRUE), 
         prepanel=function(x,y,...){
             x <- x[,drop=TRUE]
             prepanel.default.bwplot(x,y,...)},
         panel = function(x,y,mis=mes,...){
             x <- x[,drop=TRUE]
             panel.barchart(x,y,col=dots$p.bg[1],fill=dots$p.bg[1],...)
             panel.lines(type="o",mis[names(mis)==as.character(packet.number())],col=dots$l.col[1],lwd=dots$lwd[1],pch=dots$pch[1],cex=dots$cex[1],lty=dots$lty[1])},
         layout=c(1,length(unique(TechDf$GRP))),ylim=c(0,1)) 
})


                                                                                 