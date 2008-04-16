#################################################################################################
##                                                                                             ##
## Plots of volume of landings and number of fish measured by time, technical and space strata ##
##                                                                                             ##
##                                   MM 06/02/2008                                             ##
#################################################################################################


setGeneric("biasPlot", function(samObj,
                                popObj,
                                ...){
	standardGeneric("biasPlot")
})

         
                                                                     
setMethod("biasPlot", signature(samObj="csDataVal",popObj="clDataVal"), function(samObj,
                                                                                 popObj=clDataVal(),       
                                                                                 samFld="lenNum",    #or "wt","subSampWt", "nbSamp" (number of samples)
                                                                                 popFld="landWt", #or... 
                                                                                 timeStrata="quarter",  #or "year", "semester","month"
                                                                                 spaceStrata="area",    #or "rect"
                                                                                 techStrata="commCat",   #or "foCatNat" or...
                                                                                 show="all",   #or "samp", "pop"
                                                                                 ...){


CL <- FALSE
if (nrow(popObj@cl)==1) CL <- TRUE                                                                                  
TechInd <- is.null(techStrata)
TimeInd <- is.null(timeStrata)
SpaceInd <- is.null(spaceStrata)


if (all(c(TechInd,TimeInd,SpaceInd))) stop("No stratification!!")
#Tests on stratification parameters
if (!TechInd) {
  if (!techStrata%in%c("commCat","foCatNat","foCatEu5","foCatEu6")) stop("Wrong 'techStrata' parameter!")}
if (!TimeInd) {
  if (!timeStrata%in%c("year","semester","quarter","month")) stop("Wrong 'timeStrata' parameter!")}           
if (!SpaceInd) {
  if (!spaceStrata%in%c("area","rect")) stop("Wrong 'spaceStrata' parameter!")}           


#If timeStrata="semester", "quarter" or "month", field must be put in HH
HHtime <- as.numeric(as.character(samObj@hh$year))
HHmonth <- as.numeric(sapply(samObj@hh$date,function(x) strsplit(as.character(x),"-")[[1]][2]))
if (!TimeInd) {
  if (timeStrata=="month") HHtime <- HHmonth 
  if (timeStrata=="quarter") HHtime <- ceiling(HHmonth/3)
  if (timeStrata=="semester") HHtime <- ceiling(HHmonth/6)}
  
  
HHtime <- factor(HHtime,exclude=NULL)  #sorting time occurrences
if (!TimeInd) 
  eval(parse('',text=paste("samObj@hh$",timeStrata,"<- as.numeric(as.character(HHtime))",sep=""))) 
#insert 'Stratas' from HH to SL or HL 
if (samFld=="lenNum") 
  tab <- merge(samObj@hl,samObj@hh,all.x=TRUE) 
else 
  tab <- merge(samObj@sl,samObj@hh,all.x=TRUE)  

mod <- c("\"Time\"","\"Space\"","\"Technical\"")[c(!TimeInd,!SpaceInd,!TechInd)]

tabPop <- cl(popObj)
tabPop$semester <- ceiling(tabPop$quarter/2) 

if (samFld=="nbSamp"){
  fun <- "length"
  Var <- "wt"
} else {
  fun <- "sum,na.rm=TRUE"
  Var <- samFld}

index <- c(1:sum(!TimeInd,!TechInd,!SpaceInd))
eval(parse('',text=paste("Lev",index," <- unique(c(tab$",c(timeStrata,spaceStrata,techStrata),",tabPop$",                             #levels
                         c(timeStrata,spaceStrata,techStrata),"))",sep="",collapse=";")))
eval(parse('',text=paste("Mes",index," <- tapply(tab$",Var,",list(factor(tab$",c(timeStrata,spaceStrata,techStrata),                 #sampled value
                         ",levels=sort(Lev",index,"))),",fun,")",sep="",collapse=";")))
if (!CL) {eval(parse('',text=paste("Pop",index," <- tapply(tabPop$",popFld,",list(factor(tabPop$",c(timeStrata,spaceStrata,techStrata), #population value
                                   ",levels=sort(Lev",index,"))),sum,na.rm=TRUE)",sep="",collapse=";")))}



eval(parse('',text=paste("Mes",index,"[is.na(Mes",index,")] <- 0",sep="",collapse=";")))                                              #NA <- 0
if (!CL) eval(parse('',text=paste("Pop",index,"[is.na(Pop",index,")] <- 0",sep="",collapse=";")))

eval(parse('',text=paste("Mes",index," <- Mes",index,"/sum(Mes",index,",na.rm=TRUE)",sep="",collapse=";")))                           #Rate calculation
if (!CL) eval(parse('',text=paste("Pop",index," <- Pop",index,"/sum(Pop",index,",na.rm=TRUE)",sep="",collapse=";")))

eval(parse('',text=paste("dimVec <- c(",paste("length(Mes",index,")",sep="",collapse=","),")",sep="")))                               #dimensions

int1 <- paste("names(Mes",index,")",sep="",collapse=",")
int2 <- paste("Mes",index,sep="",collapse=",")
if (CL) 
  int3 <- "0" 
else 
  int3 <- paste("Pop",index,sep="",collapse=",")
int4 <- paste(mod,sep="",collapse=",") 

eval(parse('',text=paste("TechDf <- data.frame(str=c(",int1,"),mes=c(",int2,"),pop=c(",int3,"),GRP=rep(c(",int4,"),dimVec))",sep="")))  #data.frame
TechDf$str <- factor(TechDf$str,levels=as.character(TechDf$str))

TechDf$grp <- factor(TechDf$GRP,labels=1:length(unique(TechDf$GRP)))
TechDf$GRPplus <- factor(TechDf$GRP,levels=c("Space","Technical","Time"),labels=c(paste("Space =",spaceStrata),paste("Technical =",techStrata),paste("Time =",timeStrata)))
TAB <- TechDf[order(TechDf$GRP),c("str","mes","pop","GRPplus")]
rownames(TAB) <- 1:nrow(TAB)
if (show=="samp") 
  TechDf$pop <- 0 
if (show=="pop") 
  TechDf$mes <- NA

mes <- TechDf$mes
names(mes) <- TechDf$grp

#-----> graphical display 
data(GraphsPar)                                                                                                                           
dots <- list(...) 
if (is.null(dots$p.col)) 
  dots$p.col <- "black"
if (is.null(dots$p.bg)) 
  dots$p.bg <- "lightblue"
if (is.null(dots$l.col)) 
  dots$l.col <- "red"
if (is.null(dots$lwd)) 
  dots$lwd <- 2 
if (is.null(dots$pch)) 
  dots$pch <- 20
if (is.null(dots$rot)) 
  dots$rot <- 0  

sapply(names(GP),function(x) 
                  if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) 
                    eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) 
  dots$xlab <- "" 
if (is.null(dots$ylab)) 
  dots$ylab <- "Frequency" 
if (is.null(dots$main)) 
  dots$main <- paste("Relative Rates of \"",samFld,"\"",paste(c(" and \"",popFld,"\" variables"),collapse="")[!CL],"\nby ",paste(unique(TechDf$GRP),collapse=", ")," Strata",sep="")
g <- c(!CL,TRUE)
print(barchart(pop~str|GRPplus,data=TechDf,scales=list(x=list(relation="free",rot=dots$rot[1]),font=dots$font.axis),main=list(dots$main,font=dots$font.main),
               xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),par.strip.text=list(font=dots$font.lab),
               key =list(lines=list(pch=c(15,1)[g],type=c("p","l")[g],col=c(dots$p.bg[1],dots$l.col[1])[g],lwd=c(2,dots$lwd[1])[g],cex=c(1.2,dots$cex[1])[g],lty=c(1,dots$lty[1])[g]),
                         text=list(c(popFld,samFld)[g]),font=dots$font.lab,space="right",columns=1,border=TRUE), 
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





#Same thing with CE table

                                                                   
setMethod("biasPlot", signature(samObj="csDataVal",popObj="ceDataVal"), function(samObj,
                                                                                 popObj=ceDataVal(),        
                                                                                 samFld="lenNum",    #or "wt","subSampWt" , "nbSamp" (number of samples)
                                                                                 popFld="trpNum", #or ...
                                                                                 timeStrata="quarter",  #or "year", "semester","month"
                                                                                 spaceStrata="area",    #or "rect"
                                                                                 techStrata="foCatEu5",   #or "foCatNat", "foCatEu6"
                                                                                 show="all",   #or "samp", "pop"
                                                                                 ...){

CE <- FALSE
if (nrow(popObj@ce)==1) CE <- TRUE                                                                                  
TechInd <- is.null(techStrata) ; TimeInd <- is.null(timeStrata) ; SpaceInd <- is.null(spaceStrata)
if (all(c(TechInd,TimeInd,SpaceInd))) stop("No stratification!!")

#Tests on stratification parameters
if (!TechInd) {
  if (!techStrata%in%c("foCatNat","foCatEu5","foCatEu6")) stop("Wrong 'techStrata' parameter!")}
if (!TimeInd) {
  if (!timeStrata%in%c("year","semester","quarter","month")) stop("Wrong 'timeStrata' parameter!")}           
if (!SpaceInd) {
  if (!spaceStrata%in%c("area","rect")) stop("Wrong 'spaceStrata' parameter!")}           

#If timeStrata="semester", "quarter" or "month", field must be put in HH
HHtime <- as.numeric(as.character(samObj@hh$year))
HHmonth <- as.numeric(sapply(samObj@hh$date,function(x) strsplit(as.character(x),"-")[[1]][2]))
if (!TimeInd) {
  if (timeStrata=="month") HHtime <- HHmonth }
if (!TimeInd) {
  if (timeStrata=="quarter") HHtime <- ceiling(HHmonth/3) }
if (!TimeInd) {
  if (timeStrata=="semester") HHtime <- ceiling(HHmonth/6) }
  
HHtime <- factor(HHtime,exclude=NULL)  #sorting time occurrences
if (!TimeInd) eval(parse('',text=paste("samObj@hh$",timeStrata,"<- as.numeric(as.character(HHtime))",sep=""))) 
#insert 'Stratas' from HH to SL or HL 
if (samFld=="lenNum") 
  tab <- merge(samObj@hl,samObj@hh,all.x=TRUE) 
else 
  tab <- merge(samObj@sl,samObj@hh,all.x=TRUE)  

mod <- c("\"Time\"","\"Space\"","\"Technical\"")[c(!TimeInd,!SpaceInd,!TechInd)]

tabPop <- popObj@ce
tabPop$semester <- ceiling(tabPop$quarter/2) 

if (samFld=="nbSamp"){
  fun <- "length"
  Var <- "wt"
} else {
  fun <- "sum,na.rm=TRUE"
  Var <- samFld}


index <- c(1:sum(!TimeInd,!TechInd,!SpaceInd))
eval(parse('',text=paste("Lev",index," <- unique(c(tab$",c(timeStrata,spaceStrata,techStrata),",tabPop$",
                         c(timeStrata,spaceStrata,techStrata),"))",sep="",collapse=";")))
eval(parse('',text=paste("Mes",index," <- tapply(tab$",Var,",list(factor(tab$",c(timeStrata,spaceStrata,techStrata),
                         ",levels=sort(Lev",index,"))),",fun,")",sep="",collapse=";")))
if (!CE) {eval(parse('',text=paste("Pop",index," <- tapply(tabPop$",popFld,",list(factor(tabPop$",c(timeStrata,spaceStrata,techStrata),
                                   ",levels=sort(Lev",index,"))),sum,na.rm=TRUE)",sep="",collapse=";")))}



eval(parse('',text=paste("Mes",index,"[is.na(Mes",index,")] <- 0",sep="",collapse=";"))) 
if (!CE) eval(parse('',text=paste("Pop",index,"[is.na(Pop",index,")] <- 0",sep="",collapse=";")))

eval(parse('',text=paste("Mes",index," <- Mes",index,"/sum(Mes",index,",na.rm=TRUE)",sep="",collapse=";")))
if (!CE) eval(parse('',text=paste("Pop",index," <- Pop",index,"/sum(Pop",index,",na.rm=TRUE)",sep="",collapse=";")))

eval(parse('',text=paste("dimVec <- c(",paste("length(Mes",index,")",sep="",collapse=","),")",sep="")))

int1 <- paste("names(Mes",index,")",sep="",collapse=",")
int2 <- paste("Mes",index,sep="",collapse=",")
if (CE) 
  int3 <- "0" 
else 
  int3 <- paste("Pop",index,sep="",collapse=",")
int4 <- paste(mod,sep="",collapse=",") 

eval(parse('',text=paste("TechDf <- data.frame(str=c(",int1,"),mes=c(",int2,"),pop=c(",int3,"),GRP=rep(c(",int4,"),dimVec))",sep="")))
TechDf$str <- factor(TechDf$str,levels=as.character(TechDf$str))

TechDf$grp <- factor(TechDf$GRP,labels=1:length(unique(TechDf$GRP)))
TechDf$GRPplus <- factor(TechDf$GRP,levels=c("Space","Technical","Time"),labels=c(paste("Space =",spaceStrata),paste("Technical =",techStrata),paste("Time =",timeStrata)))
TAB <- TechDf[order(TechDf$GRP),c("str","mes","pop","GRPplus")]
rownames(TAB) <- 1:nrow(TAB)
if (show=="samp") 
  TechDf$pop <- 0
if (show=="pop") 
  TechDf$mes <- NA

mes <- TechDf$mes
names(mes) <- TechDf$grp

#graphical display 
data(GraphsPar)                                                                                                                           
dots <- list(...) 
if (is.null(dots$p.col)) 
  dots$p.col <- "black"
if (is.null(dots$p.bg)) 
  dots$p.bg <- "lightblue"
if (is.null(dots$l.col)) 
  dots$l.col <- "red"
if (is.null(dots$lwd)) 
  dots$lwd <- 2
if (is.null(dots$pch)) 
  dots$pch <- 20
if (is.null(dots$rot)) 
  dots$rot <- 0 

sapply(names(GP),function(x) 
                  if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) 
                    eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) 
  dots$xlab <- ""
if (is.null(dots$ylab)) 
  dots$ylab <- "Frequency" 
if (is.null(dots$main)) 
  dots$main <- paste("Relative Rates of \"",samFld,"\"",paste(c(" and \"",popFld,"\" variables"),collapse="")[!CE],"\nby ",paste(unique(TechDf$GRP),collapse=", ")," Strata",sep="")
g <- c(!CE,TRUE)
print(barchart(pop~str|GRPplus,data=TechDf,scales=list(x=list(relation="free",rot=dots$rot[1]),font=dots$font.axis),main=list(dots$main,font=dots$font.main),
               xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),par.strip.text=list(font=dots$font.lab),
               key =list(lines=list(pch=c(15,1)[g],type=c("p","l")[g],col=c(dots$p.bg[1],dots$l.col[1])[g],lwd=c(2,dots$lwd[1])[g],cex=c(1.2,dots$cex[1])[g],lty=c(1,dots$lty[1])[g]),
                         text=list(c(popFld,samFld)[g]),font=dots$font.lab,space="right",columns=1,border=TRUE), 
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

           




###################################
###################################
##  Consolidated data structure  ##
###################################
###################################



                                                                     
setMethod("biasPlot", signature(samObj="csDataCons",popObj="clDataCons"), function(samObj,
                                                                                   popObj,       
                                                                                   samFld="lenNum",    #or "wt","subSampWt" ,"nbSamp"
                                                                                   popFld="landWt", #or... 
                                                                                   show="all",   #or "samp", "pop"
                                                                                   ...){

if (samFld=="lenNum") 
  tab <- hl(samObj) 
else 
  tab <- sl(samObj)

tabPop <- cl(popObj)

#test
if (all(is.na(tab$time)) | all(is.na(tabPop$time))) {
  timeStrata <- NULL 
  TimeInd <- TRUE
} else {
  timeStrata <- "time" 
  TimeInd <- FALSE}
  
if (all(is.na(tab$space)) | all(is.na(tabPop$space))) {
  spaceStrata <- NULL 
  SpaceInd <- TRUE
} else {
  spaceStrata <- "space"
  SpaceInd <- FALSE}
  
if (all(is.na(tab$technical)) | all(is.na(tabPop$technical))) {
  techStrata <- NULL
  TechInd <- TRUE
} else {
  techStrata <- "technical" 
  TechInd <- FALSE}
                                                                                 
if (all(c(TechInd,TimeInd,SpaceInd))) stop("No stratification!!")

mod <- c("\"Time\"","\"Space\"","\"Technical\"")[c(!TimeInd,!SpaceInd,!TechInd)]


if (samFld=="nbSamp"){
  fun <- "length"
  Var <- "wt"
} else {
  fun <- "sum,na.rm=TRUE"
  Var <- samFld}




index <- c(1:sum(!TimeInd,!TechInd,!SpaceInd))
eval(parse('',text=paste("Lev",index," <- unique(c(as.character(tab$",c(timeStrata,spaceStrata,techStrata),
                         "),as.character(tabPop$",c(timeStrata,spaceStrata,techStrata),")))",sep="",collapse=";")))
eval(parse('',text=paste("Mes",index," <- tapply(tab$",Var,",list(factor(tab$",c(timeStrata,spaceStrata,techStrata),
                         ",levels=sort(Lev",index,")[sort(Lev",index,")!=\"NA\"])),",fun,")",sep="",collapse=";")))
eval(parse('',text=paste("Pop",index," <- tapply(tabPop$",popFld,",list(factor(tabPop$",c(timeStrata,spaceStrata,techStrata),",levels=sort(Lev",index,
                         ")[sort(Lev",index,")!=\"NA\"])),sum,na.rm=TRUE)",sep="",collapse=";")))



eval(parse('',text=paste("Mes",index,"[is.na(Mes",index,")] <- 0",sep="",collapse=";"))) 
eval(parse('',text=paste("Pop",index,"[is.na(Pop",index,")] <- 0",sep="",collapse=";")))

eval(parse('',text=paste("Mes",index," <- Mes",index,"/sum(Mes",index,",na.rm=TRUE)",sep="",collapse=";")))
eval(parse('',text=paste("Pop",index," <- Pop",index,"/sum(Pop",index,",na.rm=TRUE)",sep="",collapse=";")))

eval(parse('',text=paste("dimVec <- c(",paste("length(Mes",index,")",sep="",collapse=","),")",sep="")))

int1 <- paste("names(Mes",index,")",sep="",collapse=",")
int2 <- paste("Mes",index,sep="",collapse=",")
int3 <- paste("Pop",index,sep="",collapse=",")
int4 <- paste(mod,sep="",collapse=",") 

eval(parse('',text=paste("TechDf <- data.frame(str=c(",int1,"),mes=c(",int2,"),pop=c(",int3,"),GRP=rep(c(",int4,"),dimVec))",sep="")))
TechDf$str <- factor(TechDf$str,levels=as.character(TechDf$str))

TechDf$grp <- factor(TechDf$GRP,labels=1:length(unique(TechDf$GRP)))
TechDf$GRPplus <- factor(TechDf$GRP,levels=c("Space","Technical","Time"),labels=c(paste("Space =",spaceStrata),paste("Technical =",techStrata),paste("Time =",timeStrata)))
TAB <- TechDf[order(TechDf$GRP),c("str","mes","pop","GRPplus")]
rownames(TAB) <- 1:nrow(TAB)

if (show=="samp") 
  TechDf$pop <- 0 
if (show=="pop") 
  TechDf$mes <- NA

mes <- TechDf$mes
names(mes) <- TechDf$grp

#graph display 
data(GraphsPar)                                                                                                                           
dots <- list(...) 
if (is.null(dots$p.col)) 
  dots$p.col <- "black"
if (is.null(dots$p.bg)) 
  dots$p.bg <- "lightblue"
if (is.null(dots$l.col)) 
  dots$l.col <- "red"
if (is.null(dots$lwd)) 
  dots$lwd <- 2 
if (is.null(dots$pch)) 
  dots$pch <- 20
if (is.null(dots$rot)) 
  dots$rot <- 0  

sapply(names(GP),function(x) 
                  if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) 
                    eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) 
  dots$xlab <- "" 
if (is.null(dots$ylab)) 
  dots$ylab <- "Frequency" 
if (is.null(dots$main)) 
  dots$main <- paste("Relative Rates of \"",samFld,"\"",paste(c(" and \"",popFld,"\" variables"),collapse=""),"\nby ",paste(unique(TechDf$GRP),collapse=", ")," Strata",sep="")

print(barchart(pop~str|GRPplus,data=TechDf,scales=list(x=list(relation="free",rot=dots$rot[1]),font=dots$font.axis),main=list(dots$main,font=dots$font.main),
               xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),par.strip.text=list(font=dots$font.lab),
               key =list(lines=list(pch=c(15,1),type=c("p","l"),col=c(dots$p.bg[1],dots$l.col[1]),lwd=c(2,dots$lwd[1]),cex=c(1.2,dots$cex[1]),lty=c(1,dots$lty[1])),
                         text=list(c(popFld,samFld)),font=dots$font.lab,space="right",columns=1,border=TRUE), 
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



#the same with effort data

                                                                     
setMethod("biasPlot", signature(samObj="csDataCons",popObj="ceDataCons"), function(samObj,
                                                                                   popObj,       
                                                                                   samFld="lenNum",    #or "wt","subSampWt" , "nbSamp" (number of samples)
                                                                                   popFld="trpNum", #or... 
                                                                                   show="all",   #or "samp", "pop"
                                                                                   ...){

if (samFld=="lenNum") 
  tab <- hl(samObj) 
else 
  tab <- sl(samObj)

tabPop <- ce(popObj)

#test
if (all(is.na(tab$time)) | all(is.na(tabPop$time))) {
  timeStrata <- NULL
  TimeInd <- TRUE 
} else {
  timeStrata <- "time" 
  TimeInd <- FALSE}
  
if (all(is.na(tab$space)) | all(is.na(tabPop$space))) {
  spaceStrata <- NULL 
  SpaceInd <- TRUE 
} else {
  spaceStrata <- "space"
  SpaceInd <- FALSE}

if (all(is.na(tab$technical)) | all(is.na(tabPop$technical))) {
  techStrata <- NULL 
  TechInd <- TRUE
} else {
  techStrata <- "technical" 
  TechInd <- FALSE}
                                                                                 
if (all(c(TechInd,TimeInd,SpaceInd))) stop("No stratification!!")

mod <- c("\"Time\"","\"Space\"","\"Technical\"")[c(!TimeInd,!SpaceInd,!TechInd)]


if (samFld=="nbSamp"){
  fun <- "length"
  Var <- "wt"
} else {
  fun <- "sum,na.rm=TRUE"
  Var <- samFld}


index <- c(1:sum(!TimeInd,!TechInd,!SpaceInd))
eval(parse('',text=paste("Lev",index," <- unique(c(as.character(tab$",c(timeStrata,spaceStrata,techStrata),
                         "),as.character(tabPop$",c(timeStrata,spaceStrata,techStrata),")))",sep="",collapse=";")))
eval(parse('',text=paste("Mes",index," <- tapply(tab$",Var,",list(factor(tab$",c(timeStrata,spaceStrata,techStrata),
                         ",levels=sort(Lev",index,")[sort(Lev",index,")!=\"NA\"])),",fun,")",sep="",collapse=";")))
eval(parse('',text=paste("Pop",index," <- tapply(tabPop$",popFld,",list(factor(tabPop$",c(timeStrata,spaceStrata,techStrata),
                         ",levels=sort(Lev",index,")[sort(Lev",index,")!=\"NA\"])),sum,na.rm=TRUE)",sep="",collapse=";")))



eval(parse('',text=paste("Mes",index,"[is.na(Mes",index,")] <- 0",sep="",collapse=";"))) 
eval(parse('',text=paste("Pop",index,"[is.na(Pop",index,")] <- 0",sep="",collapse=";")))

eval(parse('',text=paste("Mes",index," <- Mes",index,"/sum(Mes",index,",na.rm=TRUE)",sep="",collapse=";")))
eval(parse('',text=paste("Pop",index," <- Pop",index,"/sum(Pop",index,",na.rm=TRUE)",sep="",collapse=";")))

eval(parse('',text=paste("dimVec <- c(",paste("length(Mes",index,")",sep="",collapse=","),")",sep="")))

int1 <- paste("names(Mes",index,")",sep="",collapse=",")
int2 <- paste("Mes",index,sep="",collapse=",")
int3 <- paste("Pop",index,sep="",collapse=",")
int4 <- paste(mod,sep="",collapse=",") 

eval(parse('',text=paste("TechDf <- data.frame(str=c(",int1,"),mes=c(",int2,"),pop=c(",int3,"),GRP=rep(c(",int4,"),dimVec))",sep="")))
TechDf$str <- factor(TechDf$str,levels=as.character(TechDf$str))

TechDf$grp <- factor(TechDf$GRP,labels=1:length(unique(TechDf$GRP)))
TechDf$GRPplus <- factor(TechDf$GRP,levels=c("Space","Technical","Time"),labels=c(paste("Space =",spaceStrata),paste("Technical =",techStrata),paste("Time =",timeStrata)))
TAB <- TechDf[order(TechDf$GRP),c("str","mes","pop","GRPplus")]
rownames(TAB) <- 1:nrow(TAB)

if (show=="samp") 
  TechDf$pop <- 0
if (show=="pop") 
  TechDf$mes <- NA

mes <- TechDf$mes ; names(mes) <- TechDf$grp

#graph display 
data(GraphsPar)                                                                                                                           
dots <- list(...) 
if (is.null(dots$p.col)) 
  dots$p.col <- "black"
if (is.null(dots$p.bg)) 
  dots$p.bg <- "lightblue"
if (is.null(dots$l.col)) 
  dots$l.col <- "red"
if (is.null(dots$lwd)) 
  dots$lwd <- 2 
if (is.null(dots$pch)) 
  dots$pch <- 20 
if (is.null(dots$rot)) 
  dots$rot <- 0  

sapply(names(GP),function(x) 
                  if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) 
                    eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) 
  dots$xlab <- ""
if (is.null(dots$ylab)) 
  dots$ylab <- "Frequency" 
if (is.null(dots$main)) 
  dots$main <- paste("Relative Rates of \"",samFld,"\"",paste(c(" and \"",popFld,"\" variables"),collapse=""),"\nby ",paste(unique(TechDf$GRP),collapse=", ")," Strata",sep="")

print(barchart(pop~str|GRPplus,data=TechDf,scales=list(x=list(relation="free",rot=dots$rot[1]),font=dots$font.axis),main=list(dots$main,font=dots$font.main),
               xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),par.strip.text=list(font=dots$font.lab),
               key =list(lines=list(pch=c(15,1),type=c("p","l"),col=c(dots$p.bg[1],dots$l.col[1]),lwd=c(2,dots$lwd[1]),cex=c(1.2,dots$cex[1]),lty=c(1,dots$lty[1])),
                         text=list(c(popFld,samFld)),font=dots$font.lab,space="right",columns=1,border=TRUE), 
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


                                                                 