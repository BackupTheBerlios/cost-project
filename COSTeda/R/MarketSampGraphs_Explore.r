#################################################################################################
##                                                                                             ##
## Plots of volume of landings and number of fish measured by time, technical and space strata ##
##                                                                                             ##
##                                   MM 06/02/2008                                             ##
#################################################################################################


setGeneric("SampComp.plot", function(object1,
                                     object2,
                                     var1="lenNum",
                                     var2,
                                     TimeStrat="quarter",
                                     TechStrat,
                                     SpaceStrat="area",
                                     show="all",
                                     ...){
	standardGeneric("SampComp.plot")
})

         
                                                                     
setMethod("SampComp.plot", signature(object1="csData",object2="clData"), function(object1,
                                                                                  object2=clData(),       
                                                                                  var1="lenNum",    #or "wt","subSampWt"
                                                                                  var2, #landWt or...
                                                                                  TimeStrat="quarter",  #or "year", "semester","month"
                                                                                  SpaceStrat="area",    #or "rect"
                                                                                  TechStrat,   #"commCat" or "foCatNat", "foCatEu5", "foCatEu6"
                                                                                  show="all",   #or "samp", "pop"
                                                                                  ...){

if (missing(TechStrat)) TechStrat <- "commCat" ; if (missing(var2)) var2 <- "landWt"
CL <- FALSE
if (nrow(object2@cl)==1) CL <- TRUE                                                                                  
TechInd <- is.null(TechStrat) ; TimeInd <- is.null(TimeStrat) ; SpaceInd <- is.null(SpaceStrat)
if (all(c(TechInd,TimeInd,SpaceInd))) stop("No stratification!!")
#Tests on stratification parameters
if (!TechInd) {if (!TechStrat%in%c("commCat","foCatNat","foCatEu5","foCatEu6")) stop("Wrong 'TechStrat' parameter!")}
if (!TimeInd) {if (!TimeStrat%in%c("year","semester","quarter","month")) stop("Wrong 'TimeStrat' parameter!")}           
if (!SpaceInd) {if (!SpaceStrat%in%c("area","rect")) stop("Wrong 'SpaceStrat' parameter!")}           

#If TimeStrat="semester", "quarter" or "month", field must be put in HH
HHtime <- as.numeric(as.character(object1@hh$year))
HHmonth <- as.numeric(sapply(object1@hh$date,function(x) strsplit(as.character(x),"-")[[1]][2]))
if (!TimeInd) {if (TimeStrat=="month") HHtime <- HHmonth }
if (!TimeInd) {if (TimeStrat=="quarter") HHtime <- floor((HHmonth-0.1)/3)+1 }
if (!TimeInd) {if (TimeStrat=="semester") HHtime <- floor((HHmonth-0.1)/6)+1 }
HHtime <- factor(HHtime,exclude=NULL)  #sorting time occurrences
if (!TimeInd) eval(parse('',text=paste("object1@hh$",TimeStrat,"<- as.numeric(as.character(HHtime))",sep=""))) 
#insert 'Stratas' from HH to SL or HL 
if (var1=="lenNum") tab <- merge(object1@hl,object1@hh,all.x=TRUE) else tab <- merge(object1@sl,object1@hh,all.x=TRUE)  

mod <- c("\"Time\"","\"Space\"","\"Technical\"")[c(!TimeInd,!SpaceInd,!TechInd)]

tabPop <- object2@cl
tabPop$semester <- floor((tabPop$quarter-0.1)/2)+1 


index <- c(1:sum(!TimeInd,!TechInd,!SpaceInd))
eval(parse('',text=paste("Lev",index," <- unique(c(tab$",c(TimeStrat,SpaceStrat,TechStrat),",tabPop$",c(TimeStrat,SpaceStrat,TechStrat),"))",sep="",collapse=";")))
eval(parse('',text=paste("Mes",index," <- tapply(tab$",var1,",list(factor(tab$",c(TimeStrat,SpaceStrat,TechStrat),",levels=sort(Lev",index,"))),sum,na.rm=TRUE)",sep="",collapse=";")))
if (!CL) {eval(parse('',text=paste("Pop",index," <- tapply(tabPop$",var2,",list(factor(tabPop$",c(TimeStrat,SpaceStrat,TechStrat),",levels=sort(Lev",index,
                                  "))),sum,na.rm=TRUE)",sep="",collapse=";")))}



eval(parse('',text=paste("Mes",index,"[is.na(Mes",index,")] <- 0",sep="",collapse=";"))) 
if (!CL) eval(parse('',text=paste("Pop",index,"[is.na(Pop",index,")] <- 0",sep="",collapse=";")))
eval(parse('',text=paste("Mes",index," <- Mes",index,"/sum(Mes",index,",na.rm=TRUE)",sep="",collapse=";")))
if (!CL) eval(parse('',text=paste("Pop",index," <- Pop",index,"/sum(Pop",index,",na.rm=TRUE)",sep="",collapse=";")))

eval(parse('',text=paste("dimVec <- c(",paste("length(Mes",index,")",sep="",collapse=","),")",sep="")))

int1 <- paste("names(Mes",index,")",sep="",collapse=",")
int2 <- paste("Mes",index,sep="",collapse=",")
if (CL) int3 <- "0" else int3 <- paste("Pop",index,sep="",collapse=",")
int4 <- paste(mod,sep="",collapse=",") 

eval(parse('',text=paste("TechDf <- data.frame(str=c(",int1,"),mes=c(",int2,"),pop=c(",int3,"),GRP=rep(c(",int4,"),dimVec))",sep="")))
TechDf$str <- factor(TechDf$str,levels=as.character(TechDf$str))

TechDf$grp <- factor(TechDf$GRP,labels=1:length(unique(TechDf$GRP)))
TechDf$GRPplus <- factor(TechDf$GRP,levels=c("Space","Technical","Time"),labels=c(paste("Space =",SpaceStrat),paste("Technical =",TechStrat),paste("Time =",TimeStrat)))
TAB <- TechDf[order(TechDf$GRP),c("str","mes","pop","GRPplus")] ; rownames(TAB) <- 1:nrow(TAB)
if (show=="samp") TechDf$pop <- 0 ; if (show=="pop") TechDf$mes <- NA

mes <- TechDf$mes ; names(mes) <- TechDf$grp

#graph display 
data(GraphsPar)                                                                                                                           
dots <- list(...) 
if (is.null(dots$p.col)) dots$p.col <- "black" ; if (is.null(dots$p.bg)) dots$p.bg <- "lightblue" ; if (is.null(dots$l.col)) dots$l.col <- "red"
if (is.null(dots$lwd)) dots$lwd <- 2 ; if (is.null(dots$pch)) dots$pch <- 20 ; if (is.null(dots$rot)) dots$rot <- 0  

sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "" ; if (is.null(dots$ylab)) dots$ylab <- "Frequency" 
if (is.null(dots$main)) dots$main <- paste("Relative Rates of \"",var1,"\"",paste(c(" and \"",var2,"\" variables"),collapse="")[!CL],"\nby ",paste(unique(TechDf$GRP),collapse=", ")," Strata",sep="")
g <- c(!CL,TRUE)
print(barchart(pop~str|GRPplus,data=TechDf,scales=list(x=list(relation="free",rot=dots$rot[1]),font=dots$font.axis),main=list(dots$main,font=dots$font.main),
         xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),par.strip.text=list(font=dots$font.lab),
         key =list(lines=list(pch=c(15,1)[g],type=c("p","l")[g],col=c(dots$p.bg[1],dots$l.col[1])[g],lwd=c(2,dots$lwd[1])[g],cex=c(1.2,dots$cex[1])[g],lty=c(1,dots$lty[1])[g]),
                   text=list(c(var2,var1)[g]),font=dots$font.lab,space="right",columns=1,border=TRUE), 
         prepanel=function(x,y,mis=mes,subscripts,...){
             x <- x[,drop=TRUE]
             prepanel.default.bwplot(x,y,...)
             },
         panel = function(x,y,mis=mes,subscripts,...){
             x <- x[,drop=TRUE]
             panel.barchart(x,y,col=dots$p.bg[1],fill=dots$p.bg[1],...)
             panel.lines(type="o",mis[names(mis)==as.character(packet.number())],col=dots$l.col[1],lwd=dots$lwd[1],pch=dots$pch[1],cex=dots$cex[1],lty=dots$lty[1])},
         layout=c(1,length(unique(TechDf$GRP))),ylim=c(0,max(c(TechDf$mes,TechDf$pop),na.rm=TRUE)*1.05)))
        
return(TAB) 
})



#Même chose mais avec les données population d'effort


                                                                     
setMethod("SampComp.plot", signature(object1="csData",object2="ceData"), function(object1,
                                                                                  object2=ceData(),        
                                                                                  var1="lenNum",    #or "wt","subSampWt"
                                                                                  var2, #"trpNum" or ...
                                                                                  TimeStrat="quarter",  #or "year", "semester","month"
                                                                                  SpaceStrat="area",    #or "rect"
                                                                                  TechStrat,   #or "foCatNat", "foCatEu6"
                                                                                  show="all",   #or "samp", "pop"
                                                                                  ...){

if (missing(TechStrat)) TechStrat <- "foCatEu5" ; if (missing(var2)) var2 <- "trpNum"
CE <- FALSE
if (nrow(object2@ce)==1) CE <- TRUE                                                                                  
TechInd <- is.null(TechStrat) ; TimeInd <- is.null(TimeStrat) ; SpaceInd <- is.null(SpaceStrat)
if (all(c(TechInd,TimeInd,SpaceInd))) stop("No stratification!!")
#Tests on stratification parameters
if (!TechInd) {if (!TechStrat%in%c("foCatNat","foCatEu5","foCatEu6")) stop("Wrong 'TechStrat' parameter!")}
if (!TimeInd) {if (!TimeStrat%in%c("year","semester","quarter","month")) stop("Wrong 'TimeStrat' parameter!")}           
if (!SpaceInd) {if (!SpaceStrat%in%c("area","rect")) stop("Wrong 'SpaceStrat' parameter!")}           

#If TimeStrat="semester", "quarter" or "month", field must be put in HH
HHtime <- as.numeric(as.character(object1@hh$year))
HHmonth <- as.numeric(sapply(object1@hh$date,function(x) strsplit(as.character(x),"-")[[1]][2]))
if (!TimeInd) {if (TimeStrat=="month") HHtime <- HHmonth }
if (!TimeInd) {if (TimeStrat=="quarter") HHtime <- floor((HHmonth-0.1)/3)+1 }
if (!TimeInd) {if (TimeStrat=="semester") HHtime <- floor((HHmonth-0.1)/6)+1 }
HHtime <- factor(HHtime,exclude=NULL)  #sorting time occurrences
if (!TimeInd) eval(parse('',text=paste("object1@hh$",TimeStrat,"<- as.numeric(as.character(HHtime))",sep=""))) 
#insert 'Stratas' from HH to SL or HL 
if (var1=="lenNum") tab <- merge(object1@hl,object1@hh,all.x=TRUE) else tab <- merge(object1@sl,object1@hh,all.x=TRUE)  

mod <- c("\"Time\"","\"Space\"","\"Technical\"")[c(!TimeInd,!SpaceInd,!TechInd)]

tabPop <- object2@ce
tabPop$semester <- floor((tabPop$quarter-0.1)/2)+1 


index <- c(1:sum(!TimeInd,!TechInd,!SpaceInd))
eval(parse('',text=paste("Lev",index," <- unique(c(tab$",c(TimeStrat,SpaceStrat,TechStrat),",tabPop$",c(TimeStrat,SpaceStrat,TechStrat),"))",sep="",collapse=";")))
eval(parse('',text=paste("Mes",index," <- tapply(tab$",var1,",list(factor(tab$",c(TimeStrat,SpaceStrat,TechStrat),",levels=sort(Lev",index,"))),sum,na.rm=TRUE)",sep="",collapse=";")))
if (!CE) {eval(parse('',text=paste("Pop",index," <- tapply(tabPop$",var2,",list(factor(tabPop$",c(TimeStrat,SpaceStrat,TechStrat),",levels=sort(Lev",index,
                                  "))),sum,na.rm=TRUE)",sep="",collapse=";")))}



eval(parse('',text=paste("Mes",index,"[is.na(Mes",index,")] <- 0",sep="",collapse=";"))) 
if (!CE) eval(parse('',text=paste("Pop",index,"[is.na(Pop",index,")] <- 0",sep="",collapse=";")))
eval(parse('',text=paste("Mes",index," <- Mes",index,"/sum(Mes",index,",na.rm=TRUE)",sep="",collapse=";")))
if (!CE) eval(parse('',text=paste("Pop",index," <- Pop",index,"/sum(Pop",index,",na.rm=TRUE)",sep="",collapse=";")))

eval(parse('',text=paste("dimVec <- c(",paste("length(Mes",index,")",sep="",collapse=","),")",sep="")))

int1 <- paste("names(Mes",index,")",sep="",collapse=",")
int2 <- paste("Mes",index,sep="",collapse=",")
if (CE) int3 <- "0" else int3 <- paste("Pop",index,sep="",collapse=",")
int4 <- paste(mod,sep="",collapse=",") 

eval(parse('',text=paste("TechDf <- data.frame(str=c(",int1,"),mes=c(",int2,"),pop=c(",int3,"),GRP=rep(c(",int4,"),dimVec))",sep="")))
TechDf$str <- factor(TechDf$str,levels=as.character(TechDf$str))

TechDf$grp <- factor(TechDf$GRP,labels=1:length(unique(TechDf$GRP)))
TechDf$GRPplus <- factor(TechDf$GRP,levels=c("Space","Technical","Time"),labels=c(paste("Space =",SpaceStrat),paste("Technical =",TechStrat),paste("Time =",TimeStrat)))
TAB <- TechDf[order(TechDf$GRP),c("str","mes","pop","GRPplus")] ; rownames(TAB) <- 1:nrow(TAB)
if (show=="samp") TechDf$pop <- 0 ; if (show=="pop") TechDf$mes <- NA

mes <- TechDf$mes ; names(mes) <- TechDf$grp

#graph display 
data(GraphsPar)                                                                                                                           
dots <- list(...) 
if (is.null(dots$p.col)) dots$p.col <- "black" ; if (is.null(dots$p.bg)) dots$p.bg <- "lightblue" ; if (is.null(dots$l.col)) dots$l.col <- "red"
if (is.null(dots$lwd)) dots$lwd <- 2 ; if (is.null(dots$pch)) dots$pch <- 20 ; if (is.null(dots$rot)) dots$rot <- 0 

sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "" ; if (is.null(dots$ylab)) dots$ylab <- "Frequency" 
if (is.null(dots$main)) dots$main <- paste("Relative Rates of \"",var1,"\"",paste(c(" and \"",var2,"\" variables"),collapse="")[!CE],"\nby ",paste(unique(TechDf$GRP),collapse=", ")," Strata",sep="")
g <- c(!CE,TRUE)
print(barchart(pop~str|GRPplus,data=TechDf,scales=list(x=list(relation="free",rot=dots$rot[1]),font=dots$font.axis),main=list(dots$main,font=dots$font.main),
         xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),par.strip.text=list(font=dots$font.lab),
         key =list(lines=list(pch=c(15,1)[g],type=c("p","l")[g],col=c(dots$p.bg[1],dots$l.col[1])[g],lwd=c(2,dots$lwd[1])[g],cex=c(1.2,dots$cex[1])[g],lty=c(1,dots$lty[1])[g]),
                   text=list(c(var2,var1)[g]),font=dots$font.lab,space="right",columns=1,border=TRUE), 
         prepanel=function(x,y,mis=mes,subscripts,...){
             x <- x[,drop=TRUE]
             prepanel.default.bwplot(x,y,...)
             },
         panel = function(x,y,mis=mes,subscripts,...){
             x <- x[,drop=TRUE]
             panel.barchart(x,y,col=dots$p.bg[1],fill=dots$p.bg[1],...)
             panel.lines(type="o",mis[names(mis)==as.character(packet.number())],col=dots$l.col[1],lwd=dots$lwd[1],pch=dots$pch[1],cex=dots$cex[1],lty=dots$lty[1])},
         layout=c(1,length(unique(TechDf$GRP))),ylim=c(0,max(c(TechDf$mes,TechDf$pop),na.rm=TRUE)*1.05)))
        
return(TAB) 
})

           




                                                                 