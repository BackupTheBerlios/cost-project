######################################
##                                  ##
## Delta plots (outliers detection) ##
##                                  ##
##       MM 14/01/2008              ##
######################################

##save 'Sole' & 'GraphsPar' .RData files and change the path


setwd("C:/")
require(COSTcore)


setClass("Delta.list",representation(Dlist="list"),prototype(Dlist=list()))
setClass("Delta.length",representation(Dlength="list"),prototype(Dlength=list()))



setGeneric("Delta.cs", function(object,...){
	standardGeneric("Delta.cs")
	}
)


setMethod("Delta.cs", signature(object="csData"), function(object,species,fraction,strata1,strata2=NULL,elmts1="all",elmts2="all",show.data=FALSE,...){

Hh <- hh(object)
Hh$month <- sapply(Hh$date,function(x) as.numeric(strftime(strptime(x,"%Y-%m-%d"),"%m")))
Hh$quarter <- floor((Hh$month-0.5)/3)+1
Hh$OP <- rownames(Hh)    #FO field
Hl <- hl(object) ; Hl <- Hl[Hl$spp%in%species,] ; Hl <- Hl[Hl$catchCat%in%fraction,]

fields <- c("year","quarter","month","date","time","area","rect","foDep","foCatNat","foCatEu5","foCatEu6","gear","meshSize","selDev","meshSizeSelDev")
#Test
stratArg <- match.arg(c(strata1,strata2),fields,several.ok=TRUE) 
if (!length(stratArg)%in%c(1:2)) stop("Something wrong with strata definition!!")

strata1 <- stratArg[1] ; eval(parse('',text=paste("strata2 <- ",c("NULL","stratArg[2]")[c(is.na(stratArg[2]),!is.na(stratArg[2]))],sep="")))

Hl1 <- merge(Hl,Hh[,c("sampType","landCtry","vslFlgCtry","year","proj","trpCode","staNum",strata1,"OP")],all.x=TRUE)

Hl2 <- merge(Hl1,Hh[,c("sampType","landCtry","vslFlgCtry","year","proj","trpCode","staNum",strata2)],all.x=TRUE)

Sl <- sl(object) ; Sl <- Sl[Sl$spp%in%species,] ; Sl <- Sl[Sl$catchCat%in%fraction,]
Hl3 <- merge(Hl2,Sl[,c("sampType","landCtry","vslFlgCtry","year","proj","trpCode","staNum","spp","catchCat","landCat","commCatScl","commCat","subSampCat","wt")],all.x=TRUE)

djku <- tapply(Hl3$lenNum,list(as.numeric(Hl3$OP),as.numeric(Hl3$lenCls)),sum,na.rm=TRUE) ; djku[is.na(djku)] <- 0

wHl <- unique(Hl3[,!(names(Hl3)%in%c("lenCls","lenNum"))]) ; wku <- tapply(wHl$wt,list(wHl$OP),sum,na.rm=TRUE)

#Delta matrix
Som.djku <- apply(djku,2,sum,na.rm=TRUE)
Som.wku <- sum(wku,na.rm=TRUE)
X <- (wku%*%t(Som.djku))/Som.wku
delta <- djku-X


if (!("all"%in%elmts1)) Hl3 <- Hl3[Hl3[,strata1]%in%elmts1,]
if (!("all"%in%elmts2)) Hl3 <- Hl3[Hl3[,strata2]%in%elmts2,]

strati.OP <- unique(Hl3[,c("OP",strata1,strata2)])
tri.OP <-  eval(parse('',text=paste("order(strati.OP[,strata1]",",strati.OP[,strata2]"[!is.null(strata2)],",as.numeric(strati.OP[,\"OP\"]))",sep="")))

#sorting process
strati.OP <- strati.OP[tri.OP,]

VecDelta <- apply(delta,1,sum,na.rm=TRUE)
strati.delta <- data.frame(strati.OP,Delta=VecDelta[strati.OP$OP])
row.names(strati.delta) <- 1:dim(strati.delta)[1] 
#Delta.list object creation 
res <- new("Delta.list",Dlist=list(species=species,fraction=fraction,
                                   strata=eval(parse('',text=paste("list(",strata1,"=elmts1",paste(",",strata2,"=elmts2",sep="")[!is.null(strata2)],")",sep=""))),
                                   d.jku=djku,w.ku=wku,delta.jku=delta,Delta=strati.delta))
if (show.data) return(res) else invisible(res) 
})






setMethod("plot",signature("Delta.list"), function(x,y=NULL,selection=FALSE,show.legend="right",...){  
 
 require(lattice)
 load("GraphsPar.RData")                                                                                                                           #<<<<---- to be replaced by 'data(...)'
 object <- x@Dlist
 test <- (ncol(object$Delta)==4)
 if (test) object$Delta <- object$Delta[order(object$Delta[,2],object$Delta[,3]),] else object$Delta <- object$Delta[order(object$Delta[,2]),] 
 dots <- list(...)
 sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
 if (is.null(dots$xlab)) dots$xlab <- "Sample number" ; if (is.null(dots$ylab)) dots$ylab <- "Delta values" 
 if (is.null(dots$main)) dots$main <- paste("Delta plot / Species :",paste(object$species,collapse=", "),"   Fraction :",paste(object$fraction,collapse=", "),
                                            "\n Primary strata : ",names(object$Delta)[2]) 


object$Delta[,2] <- as.factor(object$Delta[,2]) ; if (test) {object$Delta[,3] <- ff <- as.factor(object$Delta[,3])}

  XX <- 1:nrow(object$Delta) ; YY <- object$Delta$Delta ; if (test) levels(ff) <- rep(dots$p.bg,length=length(levels(ff))) else ff <- dots$p.bg[1]
  
  delimit <- tapply(object$Delta[,2],list(object$Delta[,2]),length)
	delimit <- delimit[!is.na(delimit)]
	indLab <- cumsum(delimit)

	VecDelta <- apply(object$delta.jku,1,sum,na.rm=TRUE)
	amp <- max(VecDelta)-min(VecDelta)
  decal <- rep(c(1,-1),length=length(delimit))                                                                                    #displaying process --> maybe not necessary ---------#
  concord <- function(x,vec){                                                                                                                                                          #
  if (length(vec)==1) {return(FALSE)                                                                                                                                                   #
  } else {                                                                                                                                                                             #
  if (x==1) {if ((vec[x]+vec[x+1])<(sum(vec)/6)) {return(TRUE)} else {return(FALSE)}                                                                                                   #
  } else {                                                                                                                                                                             #
	if (x==length(vec)) {if ((vec[x-1]+vec[x])<(sum(vec)/6)) {return(TRUE)} else {return(FALSE)}                                                                                         #
	} else {                                                                                                                                                                             #
	if (!(x%in%c(1,length(vec)))) {if (((vec[x]+vec[x+1])<(sum(vec)/6))|((vec[x-1]+vec[x])<(sum(vec)/6))|((vec[x-1]+vec[x]+vec[x+1])<(sum(vec)/4))) {return(TRUE)} else {return(FALSE)}  #
	} else {return(FALSE)}}}}}                                                                                                                                                           #
	decal[!sapply(1:length(delimit),concord,delimit)] <- 0                                                                                                                               #
                                                                                                                                                                #----------------------#
print(xyplot(YY~XX,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),scales=list(font=dots$font.axis),
       key=eval(parse('',text=c("NULL",paste("list(points=list(pch=dots$pch[1],fill=as.character(levels(ff)),cex=dots$p.cex[1],lwd=dots$lwd[1],col=dots$col[1]),",
                "text=list(levels(object$Delta[,3])),title=names(object$Delta)[3],cex.title=0.8,font=dots$font.lab,space=show.legend,columns=1,border=TRUE)",sep=""))[c(!test,test)])),
       panel= function(x,y,...) {
       panel.xyplot(x,y,pch=dots$pch[1],fill=as.character(ff),cex=dots$p.cex[1],lwd=dots$lwd[1],col=dots$col[1])
       panel.abline(v=indLab[-length(indLab)]+0.5,lwd=dots$l.lwd[1],lty=dots$lty[1],col=dots$l.col[1])
       panel.abline(h=0,lwd=dots$l.lwd[1],lty=dots$lty[1],col=dots$l.col[1])
       panel.text(0.5+indLab-delimit/2,min(object$Delta$Delta)+amp*0.03+amp*0.04*decal,names(indLab),col="black",font=4)}        
       ))

 
  if (selection) {
  trellis.focus("panel",1,1)
  Reponse <- panel.identify(labels=object$Delta$OP)
	#identification process
	listOP <- object$Delta[Reponse,"OP"]
	#delta matrix of identified FO
	delta.Id <- object$delta.jku[listOP,,drop=FALSE]
new("Delta.length",Dlength=list(species=object$species,fraction=object$fraction,strata=object$strata,d.jku=object$d.jku,w.ku=object$w.ku,delta.jku=object$delta.jku,
                            Delta=object$Delta,delta.Id=delta.Id,color.code=unique(as.character(ff))))  
 }})









setMethod("plot",signature("Delta.length"), function(x,y=NULL,...){
  require(lattice)
  load("GraphsPar.RData")                                                                                                                            #<<<<----to be replaced by 'data(...)'
  object <- x@Dlength
  dots <- list(...)
 sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
 if (is.null(dots$col)) dots$col <- dots$l.col ; if (is.null(dots$lwd)) dots$lwd <- dots$l.lwd[1]   #graphical default parameters
 listOP <- row.names(object$delta.Id)
 if (is.null(dots$xlab)) dots$xlab <- "Length" ; if (is.null(dots$ylab)) dots$ylab <- "Delta values" ; if (is.null(dots$main)) dots$main <- "Delta Length Frequencies" #idem

  Som.djku <- apply(object$d.jku,2,sum,na.rm=TRUE)
  Som.wku <- sum(object$w.ku,na.rm=TRUE) 
  VecMesId <- object$d.jku[listOP,,drop=FALSE]
  SommeDelta <- round(apply(object$delta.jku,1,sum)[listOP],2)
  VecSomme <- object$w.ku[listOP]%*%t(Som.djku/Som.wku) ; if (length(listOP)==1) dimnames(VecSomme)[[1]] <- listOP
  nam <- dimnames(VecSomme)
  df <- data.frame(x1=rep(colnames(VecSomme),each=nrow(VecSomme)),y1=rep(rownames(VecSomme),ncol(VecSomme)),Ech=as.vector(VecMesId),Exp=as.vector(VecSomme))
  xyplot(Ech+Exp~x1|y1 ,data=df,type=c("h","l"),lty=rep(dots$lty,length=2),par.strip.text=list(font=dots$font.lab),
  col=rep(dots$col,length=2),lwd=rep(dots$lwd,length=2),distribute.type=TRUE,scales=list(font=dots$font.axis,x=list(rot=90)),
  main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab))
})



setGeneric("plot.Samp", function(x,...){
	standardGeneric("plot.Samp")
	}
)

setMethod("plot.Samp",signature("Delta.list"), function(x,SampNum,show.legend="right",...){
  require(lattice)
  load("GraphsPar.RData")                                                                                                                            #<<<<----to be replaced by 'data(...)'
  object <- x@Dlist
  dots <- list(...)
 sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
 if (is.null(dots$col)) dots$col <- dots$l.col ; if (is.null(dots$lwd)) dots$lwd <- dots$l.lwd[1]   #graphical default parameters
 if (is.null(dots$xlab)) dots$xlab <- "Length" ; if (is.null(dots$ylab)) dots$ylab <- "Delta values" 
 if (is.null(dots$main)) dots$main <- paste("Delta Length Frequencies  /  Sample n°",SampNum) #idem

  Som.djku <- apply(object$d.jku,2,sum,na.rm=TRUE)
  Som.wku <- sum(object$w.ku,na.rm=TRUE) 
  VecMesId <- object$d.jku[SampNum,]
  SommeDelta <- round(apply(object$delta.jku,1,sum)[SampNum],2)
  VecSomme <- object$w.ku[SampNum]%*%t(Som.djku/Som.wku)
  df <- data.frame(x1=colnames(VecSomme),Ech=as.vector(VecMesId),Exp=as.vector(VecSomme))
  xyplot(Ech+Exp~x1,data=df,type=c("h","l"),lty=rep(dots$lty,length=2),par.strip.text=list(font=dots$font.lab),col=rep(dots$col,length=2),
  lwd=rep(dots$lwd,length=2),distribute.type=TRUE,scales=list(font=dots$font.axis,x=list(rot=90)),
  key=list(lines=list(lty=rep(dots$lty,length=2),col=rep(dots$col,length=2),lwd=rep(dots$lwd,length=2)),text=list(c("Sampled","Overall")),font=dots$font.lab,space=show.legend,columns=1,border=TRUE),
  main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab))
})




###########
# Example #
###########


load("Sole.RData")                                                                                                                            #<<<<---- to be replaced by 'data(...)'
object <- sole3.cs
#only sea sampling data is kept
object@tr <- object@tr[object@tr$sampType=="S",] 
object@hh <- object@hh[object@hh$sampType=="S",] 
object@sl <- object@sl[object@sl$sampType=="S",]
object@hl <- object@hl[object@hl$sampType=="S",]


aa1 <- Delta.cs(object,"SOL","LAN","quarter")
aa2 <- Delta.cs(object,"SOL","LAN","quarter","month")
plot(aa1)
bb <- plot(aa2,selection=TRUE,show.legend="right")         #select FOs (points) and press right-click to end

plot(bb)
plot.Samp(aa2,"194")




 
 
 

 
 