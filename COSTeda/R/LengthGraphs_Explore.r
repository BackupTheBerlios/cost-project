######################################
##                                  ##
## Delta plots (outliers detection) ##
##    Length Distributions (hl)     ##
##                                  ##
##         MM 31/01/2008            ##
######################################




############################################################
# Procédure d'agrégation rapide pour remplacer 'aggregate' #
############################################################



SpeedAgreg2 <- function(X,BY,FUN,...){
FactCar <- sapply(BY,as.character)
val <- apply(FactCar,1,function(x) paste(x,collapse="::"))
valAg <- aggregate(X,list(val=val),FUN,...)
tab <- as.data.frame(matrix(unlist(strsplit(as.character(valAg$val),"::")),ncol=length(BY),byrow=TRUE))
tab.ag <- data.frame(tab,valAg[,-1])
namBY <- names(BY) ; namX <- names(X)
if (is.null(namBY)) namBY <- rep("",length(BY)) ; if (is.null(namX)) namX <- rep("",length(X))
namBY[namBY==""] <- paste("c.",1:sum(namBY==""),sep="") ; namX[namX==""] <- paste("v.",1:sum(namX==""),sep="")
names(tab.ag) <- c(namBY,namX)
return(tab.ag)}





setGeneric("UE", function(object,         #csData()
                          species,
                          ...){
	standardGeneric("UE")}
)




setMethod("UE", signature(object="csData"), function(object,
                                                     species,
                                                     ...){

hh <- object@hh ; sl <- object@sl ; hl <- object@hl

#on ne conserve que les échantillons valides
hh <- hh[hh$foVal=="V",]

procQUARTER <- function(x) {return(ceiling(as.numeric(strsplit(as.character(x),"-")[[1]][2])/3))}
procMONTH <- function(x) {return(as.numeric(strsplit(as.character(x),"-")[[1]][2]))}
hh$quarter <- sapply(hh$date,procQUARTER)
hh$month <- sapply(hh$date,procMONTH)

newsl <- sl[sl$spp%in%species,]
slhh <- merge(newsl,hh,by=c("sampType","landCtry","vslFlgCtry","year","proj","trpCode","staNum"),all=FALSE)

newhl <- hl[hl$spp%in%species,]
hlslhh <- merge(newhl,slhh,by=c("sampType","landCtry","vslFlgCtry","year","proj","trpCode","staNum","spp","catchCat","landCat","commCatScl","commCat","subSampCat"),all=FALSE)
invisible(list(spp=species,slhh=slhh,hlslhh=hlslhh))
})





############################################################################################################
############################################################################################################



setClass("DeltA",representation(Parameters="list",
                                     OutP="list"),
                      prototype(Parameters=list(),
                                OutP=list()))

setClass("DeltaID",representation(OutP="list"),
                      prototype(OutP=list()))




setGeneric("Delta", function(object,
                                   species,
                                   tempStrata=NULL, #ex: "year","quarter","month"
                                   spaceStrata=NULL, #ex: "area","rect"
                                   techStrata=NULL, #ex: "commCat","foCatEu5"
                                   indSamp=FALSE,
                                   strategy="metier",  #ou "cc"
                                   ...){
	standardGeneric("Delta")}
)




setMethod("Delta", signature(object="csData"),function(object,
                                                             species,
                                                             tempStrata=NULL, #ex: "year","quarter","month"
                                                             spaceStrata=NULL, #ex: "area","rect"
                                                             techStrata=NULL, #ex: "commCat","foCatEu5"
                                                             indSamp=FALSE,
                                                             strategy="metier", #ou "cc"
                                                             ...){

tab <- UE(object,species)
tabHL <- tab$hlslhh 
#tabHL <- tabHL[!apply(cbind(tabHL[,tempStrata],tabHL[,spaceStrata],tabHL[,techStrata]),1,function(x) any(is.na(x))),]   
tabHL$Number <- tabHL$lenNum*(tabHL$wt/tabHL$subSampWt)
#on rassemble les sous-échantillons des CC pour le calcul des variances
tb <- unique(tabHL[,c(1:13,17:18,24)]) ; tb <- tb[,-13] #on ôte les sousCC et on agrège
tabHL.Wt <- SpeedAgreg2(list(wt=tb$wt,subSampWt=tb$subSampWt),list(sampType=tb$sampType,landCtry=tb$landCtry,vslFlgCtry=tb$vslFlgCtry,year=tb$year,proj=tb$proj,
                                                                   trpCode=tb$trpCode,staNum=tb$staNum,spp=tb$spp,catchCat=tb$catchCat,landCat=tb$landCat,commCatScl=tb$commCatScl,
                                                                   commCat=tb$commCat,date=tb$date),sum,na.rm=TRUE)
#on fait de même avec les nombres aux tailles
th <- tabHL[,c(1:12,15,44,19,24)]
tabHL.Lg <- SpeedAgreg2(list(Number=th$Number),list(sampType=th$sampType,landCtry=th$landCtry,vslFlgCtry=th$vslFlgCtry,year=th$year,proj=th$proj,
                                                    trpCode=th$trpCode,staNum=th$staNum,spp=th$spp,catchCat=th$catchCat,landCat=th$landCat,commCatScl=th$commCatScl,
                                                    commCat=th$commCat,date=th$date,lenCode=th$lenCode,lenCls=th$lenCls),sum,na.rm=TRUE)
#et on rassemble le tout
un <- unique(tabHL[,c(1:12,20:(ncol(tabHL)-1))]) ; deux <- merge(un,tabHL.Wt,all.y=TRUE) ; tabHL <- merge(deux,tabHL.Lg,all.y=TRUE)
tabHL$lenCls <- as.numeric(as.character(tabHL$lenCls))
tabHL$wt <- as.numeric(as.character(tabHL$wt))
tabHL$subSampWt <- as.numeric(as.character(tabHL$subSampWt))
tabHL$Number <- as.numeric(as.character(tabHL$Number))

Ntp <- is.null(tempStrata) ; Nsp <- is.null(spaceStrata) ; Ntc <- is.null(techStrata)

if (!Ntp) {
  if (all(is.na(tabHL[,tempStrata]))) {tempStrata <- NULL ; Ntp <- TRUE}
}

if (!Nsp) {
  if (all(is.na(tabHL[,spaceStrata]))) {spaceStrata <- NULL ; Nsp <- TRUE}
}

if (!Ntc) {
  if (all(is.na(tabHL[,techStrata]))) {techStrata <- NULL ; Ntc <- TRUE}
}

if (!Ntp) tabHL[,tempStrata] <- factor(tabHL[,tempStrata],exclude=c("NA",NA))
if (!Nsp) tabHL[,spaceStrata] <- factor(tabHL[,spaceStrata],exclude=c("NA",NA))
if (!Ntc) tabHL[,techStrata] <- factor(tabHL[,techStrata],exclude=c("NA",NA))


if (strategy=="cc") {UniteInt <- apply(tabHL[,1:13],1,paste,collapse="::")                       
} else {                                                                                              
UniteInt <- apply(tabHL[,1:11],1,paste,collapse="::")}                                            
    
lev <- levels(factor(UniteInt))
tabHL$Unite <- factor(UniteInt,levels=lev,labels=1:length(lev))

#afin de retrouver l'échantillon --> on garde trpCode/staNum/spp/cc si strategy=="cc" ou trpCode/staNum/spp si strategy=="metier"
extract <- function(vec,elmts) unlist(lapply(strsplit(vec,"::"),function(x) x[elmts])) 

if (strategy=="cc") {
DFsamp <- data.frame(SampNum=(1:length(lev)),trpCode=extract(lev,6),staNum=extract(lev,7),spp=extract(lev,8),commCat=extract(lev,12))
} else {
DFsamp <- data.frame(SampNum=(1:length(lev)),trpCode=extract(lev,6),staNum=extract(lev,7),spp=extract(lev,8))
}


lenSp <- c(1,5,10,25) ; names(lenSp) <- c("mm","scm","cm","25mm")
tabHL$Length <- factor(tabHL$lenCls,levels=seq(min(tabHL$lenCls,na.rm=TRUE),max(tabHL$lenCls,na.rm=TRUE),by=lenSp[as.character(unique(tabHL$lenCode)[1])]))

eval(parse('',text=paste("nbTot_Lg <- tapply(tabHL$Number,list(tabHL$Length",","[any(c(!Ntp,!Nsp,!Ntc))],
                          paste(c("as.character(tabHL[,tempStrata])"[!Ntp],"as.character(tabHL[,spaceStrata])"[!Nsp],"as.character(tabHL[,techStrata])"[!Ntc]),collapse=",",sep=""),
                          "),sum,na.rm=TRUE)",sep="")))
TotLg <- apply(nbTot_Lg,1,sum,na.rm=TRUE)                             #nombre total sur toutes les strates croisées

if (strategy=="cc") {
eval(parse('',text=paste("WkvTot <- tapply(tabHL$wt,list(as.character(tabHL$Unite)",","[any(c(!Ntp,!Nsp,!Ntc))],
                          paste(c("as.character(tabHL[,tempStrata])"[!Ntp],"as.character(tabHL[,spaceStrata])"[!Nsp],"as.character(tabHL[,techStrata])"[!Ntc]),collapse=",",sep=""),
                          "),mean,na.rm=TRUE)",sep="")))                  
  } else {                                                                                           
  tabtabHL <- cbind(tabHL[,1:13],tabHL[,c(tempStrata,spaceStrata,techStrata,"wt","Unite")]) ; tabtabHL <- unique(tabtabHL)                       
eval(parse('',text=paste("WkvTot <- tapply(tabtabHL$wt,list(as.character(tabtabHL$Unite)",","[any(c(!Ntp,!Nsp,!Ntc))],
                          paste(c("as.character(tabtabHL[,tempStrata])"[!Ntp],"as.character(tabtabHL[,spaceStrata])"[!Nsp],"as.character(tabtabHL[,techStrata])"[!Ntc]),collapse=",",sep=""),
                          "),sum,na.rm=TRUE)",sep="")))                          
  }       

TotW <- sum(WkvTot,na.rm=TRUE)

progint <- function(x,ind,tabHL) {      #ind=1 --> delta , ind=2 --> Nk , ind=3 --> Wk , ind=4 --> DELTA
  tab <- tabHL[x,]
  Djkv <- tapply(tab$Number,list(tab$Length,as.character(tab$Unite)),sum,na.rm=TRUE) ; Djkv[is.na(Djkv)] <- 0
  
  if (strategy=="cc") {Wkv <- tapply(tab$wt,list(as.character(tab$Unite)),mean)                   
  } else {                                                                                           
  tabtab <- cbind(tab[,1:13],tab[,c("wt","Unite")]) ; tabtab <- unique(tabtab)                       
  Wkv <- tapply(tabtab$wt,list(as.character(tabtab$Unite)),sum,na.rm=TRUE)                           
  }                                                                                                  
  DjkvSum <- apply(Djkv,1,sum,na.rm=TRUE)
  WkvSum <- sum(Wkv,na.rm=TRUE)
  Nk <- length(Wkv) ; Wk <- sum(Wkv,na.rm=TRUE)
  delta <- Djkv - (DjkvSum/WkvSum)%*%t(Wkv)
  delta2 <- Djkv - (TotLg/TotW)%*%t(Wkv)
  deltaSamp <- apply(delta2,2,sum,na.rm=TRUE)
  eval(parse('',text=paste("DELTA <- data.frame(samp=names(deltaSamp),delta=deltaSamp",","[any(c(!Ntp,!Nsp,!Ntc))],
                           paste(c("tp=as.character(tab[1,tempStrata])"[!Ntp],"sp=as.character(tab[1,spaceStrata])"[!Nsp],"tc=as.character(tab[1,techStrata])"[!Ntc]),
                           collapse=",",sep=""),")",sep="")))
  deltaSq <- apply(delta^2,1,sum,na.rm=TRUE)
  ll <- list(deltaSq,Nk,Wk,DELTA)[[ind]]
  return(ll)
}


if (!indSamp) {

eval(parse('',text=paste("Delta <- tapply(1:nrow(tabHL),list(",
                          paste(c("as.character(tabHL[,tempStrata])"[!Ntp],"as.character(tabHL[,spaceStrata])"[!Nsp],"as.character(tabHL[,techStrata])"[!Ntc]),collapse=",",sep=""),
                          "),function(x) progint(x,1,tabHL))",sep="")))

DeltaMatrix <- array(0,dim=c(length(levels(tabHL$Length)),dim(Delta)),dimnames=unlist(list(list(levels(tabHL$Length)),dimnames(Delta)),recursive=FALSE))
index <- apply(expand.grid(dimnames(Delta)),1,function(x) paste("\"",x,"\"",collapse=",",sep=""))
invisible(sapply(index,function(x) eval(parse('',text=paste("if (!is.null(Delta[",index,"][[1]])) DeltaMatrix[,",index,"] <<- Delta[",index,"][[1]]",sep="")))))

eval(parse('',text=paste("NkMatrix <- tapply(1:nrow(tabHL),list(",
                          paste(c("as.character(tabHL[,tempStrata])"[!Ntp],"as.character(tabHL[,spaceStrata])"[!Nsp],"as.character(tabHL[,techStrata])"[!Ntc]),collapse=",",sep=""),
                          "),function(x) progint(x,2,tabHL))",sep="")))
eval(parse('',text=paste("WkMatrix <- tapply(1:nrow(tabHL),list(",
                          paste(c("as.character(tabHL[,tempStrata])"[!Ntp],"as.character(tabHL[,spaceStrata])"[!Nsp],"as.character(tabHL[,techStrata])"[!Ntc]),collapse=",",sep=""),
                          "),function(x) progint(x,3,tabHL))",sep="")))


invisible(new("DeltA",Parameters=list(species=species,tempStrata=tempStrata,spaceStrata=spaceStrata,techStrata=techStrata),
                  OutP=list(DeltaMatrix=DeltaMatrix,NkMatrix=NkMatrix,WkMatrix=WkMatrix)))
} else {

eval(parse('',text=paste("SampDeltaMat <- do.call(\"rbind\",tapply(1:nrow(tabHL),list(",
                          paste(c("as.character(tabHL[,tempStrata])"[!Ntp],"as.character(tabHL[,spaceStrata])"[!Nsp],"as.character(tabHL[,techStrata])"[!Ntc]),collapse=",",sep=""),
                          "),function(x) progint(x,4,tabHL)))",sep="")))

invisible(new("DeltA",Parameters=list(species=species,tempStrata=tempStrata,spaceStrata=spaceStrata,techStrata=techStrata),
                  OutP=list(SampDeltaMat=SampDeltaMat,tab=tabHL,DFsamp=DFsamp)))
}

})




setGeneric("plot.Delta", function(x,
                                  y=NULL,
                                  species,
                                  tempStrata=NULL, #ex: "year","quarter","month"
                                  spaceStrata=NULL, #ex: "area","rect"
                                  techStrata=NULL, #ex: "commCat","foCatEu5"
                                  elmts=list(tp="all",sp="all",tc="all"),
                                  strategy="metier", 
                                  strat1,strat2="NULL",# pour choisir la stratification de représentation (2 au maximum) entre tempStrata, spaceStrata et techStrata
                                  selection=FALSE,
                                  show.legend="right",
                                  ...){
	standardGeneric("plot.Delta")}
)




setMethod("plot.Delta",signature("csData"), function(x,
                                               y=NULL,
                                               species,
                                               tempStrata=NULL, #ex: "year","quarter","month"
                                               spaceStrata=NULL, #ex: "area","rect"
                                               techStrata=NULL, #ex: "commCat","foCatEu5"
                                               elmts=list(tp="all",sp="all",tc="all"),   #restriction à certaines modalités de certaines strates
                                               strategy="metier", 
                                               strat1,strat2="NULL",# pour choisir la stratification de représentation (2 au maximum) entre tempStrata, spaceStrata et techStrata
                                               selection=FALSE,
                                               show.legend="right",
                                               shift=FALSE,         #décalage des modalités selon la dimension
                                               ...){  
 


INDSAMP <- Delta(x,species=species,tempStrata=tempStrata,spaceStrata=spaceStrata,techStrata=techStrata,indSamp=TRUE,strategy=strategy)@OutP
#INDLENG <- Delta(x,species=species,tempStrata=tempStrata,spaceStrata=spaceStrata,techStrata=techStrata,indSamp=FALSE,strategy=strategy)@OutP

object <- INDSAMP$SampDeltaMat
#on opère la sélection formulée par le paramètre 'elmts'
invisible(sapply(names(object)[3:ncol(object)],function(x) {elm <- unlist(elmts[x]) ; if (!"all"%in%elm) object <<- object[as.character(object[,x])%in%elm,]})) 


#on renomme                            
if (ncol(object)>2) names(object)[3:ncol(object)] <- c(tempStrata,spaceStrata,techStrata)   #


data(GraphsPar)                                                                                                                          
dots <- list(...)
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Sample number" ; if (is.null(dots$ylab)) dots$ylab <- "Delta values" 
if (is.null(dots$main)) dots$main <- paste("Delta plot / Species :",paste(species,collapse=", "),
                                            "\n Primary strata : ",eval(parse('',text=strat1)),sep="") 


FStr <- eval(parse('',text=strat1)) ; SStr <- eval(parse('',text=strat2))
if (!is.null(SStr))  object2 <- object[order(object[,FStr],object[,SStr]),] else object2 <- object[order(object[,FStr]),]
object2[,FStr] <- factor(object2[,FStr]) ; if (!is.null(SStr)) {object2[,SStr] <- ff <- factor(object2[,SStr])}

XX <- 1:nrow(object2) ; YY <- object2$delta
if (!is.null(SStr)) levels(ff) <- rep(dots$p.bg,length=length(levels(ff))) else ff <- dots$p.bg[1]

delimit <- tapply(object2[,FStr],list(object2[,FStr]),length)
delimit <- delimit[!is.na(delimit)]
indLab <- cumsum(delimit)


amp <- max(object2$delta)-min(object2$delta)
if (shift) decal <- rep(c(1,-1),length=length(delimit)) else decal <- 0                                                                                  
                                                                                                                    
print(xyplot(YY~XX,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),scales=list(font=dots$font.axis),
                key=eval(parse('',text=c("NULL",paste("list(points=list(pch=dots$pch[1],fill=as.character(levels(ff)),cex=dots$p.cex[1],lwd=dots$lwd[1],col=dots$col[1]),",
                "text=list(levels(object2[,SStr])),title=SStr,cex.title=0.8,font=dots$font.lab,space=show.legend,columns=1,border=TRUE)",sep=""))[c(is.null(SStr),!is.null(SStr))])),
                panel= function(x,y,...) {
                panel.xyplot(x,y,pch=dots$pch[1],fill=as.character(ff),cex=dots$p.cex[1],lwd=dots$lwd[1],col=dots$col[1])
                panel.abline(v=indLab[-length(indLab)]+0.5,lwd=dots$l.lwd[1],lty=dots$lty[1],col=dots$l.col[1])
                panel.abline(h=0,lwd=dots$l.lwd[1],lty=dots$lty[1],col=dots$l.col[1])
                panel.text(0.5+indLab-delimit/2,min(object2$delta)+amp*0.03+amp*0.04*decal,names(indLab),col="black",font=4)}        
))

 

if (selection) {
  trellis.focus("panel",1,1)
  Reponse <- panel.identify(labels=object2$samp)
	#identification process
	listOP <- object2$samp[Reponse]
	#valeurs qui serviront au graphe des échantillons identifiés
	
  new("DeltaID",OutP=list(species=species,sampId=INDSAMP$DFsamp[INDSAMP$DFsamp$SampNum%in%listOP,,drop=FALSE],tabId=INDSAMP$tab[INDSAMP$tab$Unite%in%listOP,,drop=FALSE],
                       tab=INDSAMP$tab)) 
 }

})






setMethod("plot",signature("DeltaID"), function(x,
                                                y=NULL,
                                                show.legend="right",...){

tabID <- x@OutP$tabId ; TAB <- x@OutP$tab
tabID$Unite <- as.character(tabID$Unite) 
VecMesId <- tapply(tabID$Number,list(tabID$Unite,tabID$Length),sum,na.rm=TRUE,drop=FALSE)
Som.djku <- apply(tapply(TAB$Number,list(TAB$Unite,TAB$Length),sum,na.rm=TRUE,drop=FALSE),2,sum,na.rm=TRUE)
tabUn <- unique(tabID[,c("Unite","wt")])
tabUnPop <- unique(TAB[,c("Unite","wt")])
Som.wku <- sum(tabUnPop$wt,na.rm=TRUE)
w.ku <- tapply(tabUn$wt,list(tabUn$Unite),sum,na.rm=TRUE)
VecSomme <- w.ku%*%t(Som.djku/Som.wku) 

data(GraphsPar)                                                                                                                        
dots <- list(...)
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$col)) dots$col <- dots$l.col ; if (is.null(dots$lwd)) dots$lwd <- dots$l.lwd[1]  
if (is.null(dots$xlab)) dots$xlab <- "Length" ; if (is.null(dots$ylab)) dots$ylab <- "Numbers" ; if (is.null(dots$main)) dots$main <- "Delta Length Frequencies" 


df <- data.frame(x1=as.numeric(rep(colnames(VecSomme),each=nrow(VecSomme))),y1=rep(rownames(VecSomme),ncol(VecSomme)),Ech=as.vector(VecMesId),Exp=as.vector(VecSomme))     #correction sur la classe de x1

xyplot(Ech+Exp~x1|y1 ,data=df,type=c("h","l"),lty=rep(dots$lty,length=2),par.strip.text=list(font=dots$font.lab),
       col=rep(dots$col,length=2),lwd=rep(dots$lwd,length=2),distribute.type=TRUE,scales=list(font=dots$font.axis,x=list(rot=90)),
       key=list(lines=list(lty=rep(dots$lty,length=2),col=rep(dots$col,length=2),lwd=rep(dots$lwd,length=2)),text=list(c("Sampled","Overall")),
       font=dots$font.lab,space=show.legend,columns=1,border=TRUE),
       main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab))
})






setGeneric("plot.Samp", function(x,
                                 SampNum,
                                 show.legend="right",...){
	standardGeneric("plot.Samp")
})




setMethod("plot.Samp",signature("DeltaID"), function(x,
                                                     SampNum,
                                                     show.legend="right",...){



tabID <- x@OutP$tabId ; tabID <- tabID[tabID$Unite==SampNum,] ; TAB <- x@OutP$tab
tabID$Unite <- as.character(tabID$Unite) 
VecMesId <- tapply(tabID$Number,list(tabID$Unite,tabID$Length),sum,na.rm=TRUE,drop=FALSE)
Som.djku <- apply(tapply(TAB$Number,list(TAB$Unite,TAB$Length),sum,na.rm=TRUE,drop=FALSE),2,sum,na.rm=TRUE)
tabUn <- unique(tabID[,c("Unite","wt")])
tabUnPop <- unique(TAB[,c("Unite","wt")])
Som.wku <- sum(tabUnPop$wt,na.rm=TRUE)
w.ku <- tapply(tabUn$wt,list(tabUn$Unite),sum,na.rm=TRUE)
VecSomme <- w.ku%*%t(Som.djku/Som.wku)


data(GraphsPar)                                                                                                                          
dots <- list(...)
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$col)) dots$col <- dots$l.col ; if (is.null(dots$lwd)) dots$lwd <- dots$l.lwd[1]  
if (is.null(dots$xlab)) dots$xlab <- "Length" ; if (is.null(dots$ylab)) dots$ylab <- "Numbers" 
if (is.null(dots$main)) dots$main <- paste("Delta Length Frequencies  /  Sample",SampNum) 


df <- data.frame(x1=as.numeric(colnames(VecSomme)),Ech=as.vector(VecMesId),Exp=as.vector(VecSomme))          #correction sur la classe de x1

xyplot(Ech+Exp~x1,data=df,type=c("h","l"),lty=rep(dots$lty,length=2),par.strip.text=list(font=dots$font.lab),col=rep(dots$col,length=2),
       lwd=rep(dots$lwd,length=2),distribute.type=TRUE,scales=list(font=dots$font.axis,x=list(rot=90)),
       key=list(lines=list(lty=rep(dots$lty,length=2),col=rep(dots$col,length=2),lwd=rep(dots$lwd,length=2)),text=list(c("Sampled","Overall")),
       font=dots$font.lab,space=show.legend,columns=1,border=TRUE),main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),
       ylab=list(dots$ylab,font=dots$font.lab))
})











setGeneric("plot.LD", function(x,
                               species,
                               fraction="LAN",
                               trpCode,
                               staNum="all",...){
	standardGeneric("plot.LD")
})




setMethod("plot.LD",signature("csData"), function(x,
                                                  species,
                                                  fraction="LAN",
                                                  trpCode,
                                                  staNum="all",...){      #nouvelle option à rajouter "allSum"
                                                  
if (length(species)!=1) stop("Only one species!!") ; if (length(trpCode)!=1) stop("Only one trip!!")
trpCode <- as.character(trpCode) ; staNum <- as.character(staNum) ; if ("all"%in%staNum) staNum <- "all" 
data(GraphsPar)                                                                                                                       
object <- x@hl ; lgthCode <- as.character(x@sl[(x@sl$trpCode%in%trpCode)&(x@sl$spp%in%species),"lenCode"][1])
stepp <- c(1,5,10,25) ; names(stepp) <- c("mm","scm","cm","25mm") ; ste <- stepp[lgthCode]

dots <- list(...) ; if (is.null(dots$p.col)) dots$p.col <- "black" ; if (is.null(dots$p.bg)) dots$p.bg <- "lightblue" ; if (is.null(dots$cex.axis)) dots$cex.axis <- 0.8
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Length" ; if (is.null(dots$ylab)) dots$ylab <- "Number" 
if (is.null(dots$main)) dots$main <- paste("Length Distributions by samples for trip",trpCode) 

df <- object[(object$trpCode%in%trpCode)&(object$catchCat%in%fraction)&(object$spp%in%species),]
if (nrow(df)==0) stop("No data for specified trip code and fraction!!")
if ("all"%in%staNum) {staNum <- unique(as.character(df$staNum))
} else {                                                          #ajout
if ("allSum"%in%staNum) df$staNum <- staNum <- "all" }            #<----

df <- df[df$staNum%in%staNum,] ; df$staNum <- factor(df$staNum)
if (nrow(df)==0) stop("No data for specified station number!!")
#empty length classes are considered
df$lenCls <- factor(df$lenCls,levels=seq(min(df$lenCls),max(df$lenCls),by=ste))  
LD <- tapply(df$lenNum,list(staNum=df$staNum,lenCls=df$lenCls),sum,na.rm=TRUE)
LD[is.na(LD)] <- 0  ; ll <- dimnames(LD)
DF <- data.frame(staNum=rep(ll$staNum,ncol(LD)),lenCls=rep(ll$lenCls,each=nrow(LD)),val=as.numeric(LD))
if (!DF$staNum[1]=="all") DF$staNum <- factor(DF$staNum,levels=as.character(sort(as.numeric(levels(DF$staNum)))))

barchart(val~lenCls|staNum,data=DF,ylim=c(0,max(DF$val)*1.05),scales=list(x=list(rot=dots$rot,cex=dots$cex.axis),
         font=dots$font.axis),main=list(dots$main,font=dots$font.main),
         xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),
         par.strip.text=list(font=dots$font.lab),col=dots$p.bg,fill=dots$p.bg)  
})
  

