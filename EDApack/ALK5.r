################################################
##                                            ##
## Age at Length Analysis (multinomial model) ##
##                                            ##
##       MM 14/01/2008                        ##
################################################

##save 'sole2' & 'GP' .RData files and change the path
setwd("C:/")

library(COSTcore)


setClass("ALmult",representation(obj.multi="list"),prototype(obj.multi=list()))


setGeneric("AL.multi", function(object,...){
	standardGeneric("AL.multi")}
)


#no options for multinomial modelisation for the moment  ----> to be continued
       

setMethod("AL.multi", signature(object="csData"), function(object,strata,elmts="all",age.plus,...){
require(nnet)
tab <- object@ca
tab$age <- as.numeric(as.character(tab$age))   
if (is.null(strata)) {tab$nostrata <- "all" ; strata <- "nostrata"}
#grouping ages
tab$age[tab$age>=age.plus] <- age.plus
newtab <- tab[,c("age","lenCls",strata)] ; if (!("all"%in%elmts)) newtab <- newtab[as.character(newtab[,strata])%in%elmts,]
newtab[,strata] <- factor(newtab[,strata])
alk <- tapply(newtab$age,list(as.numeric(newtab$lenCls),as.numeric(newtab$age),newtab[,strata]),length)
alk[is.na(alk)] <- 0
matTest <- apply(alk,c(1,3),function(x) all(x==0)) ; test <- matrix(matTest,ncol=1)
Age <- apply(alk,2,rbind)
L1 <- dim(alk)[1] ; L2 <- dim(alk)[2] ; L3 <- dim(alk)[3]
Length <- rep(as.numeric(dimnames(alk)[[1]]),L3)
eval(parse('',text=paste(strata,"<- (rep(dimnames(alk)[[3]],each=L1))")))
eval(parse('',text=paste("Dd <- cbind.data.frame(as.data.frame(Age),data.frame(Length=Length,",strata,"=",strata,"))",sep="")))
Dd1 <- Dd[!test,] ; Age <- Age[!test,]
ntest <- strata=="nostrata"
eval(parse('',text=paste("Mm <- multinom(Age~Length",c(paste("*",strata,"-1",sep=""),"")[c(!ntest,ntest)],",Hess=T,model=T,data=Dd1)",sep="")))
return(new("ALmult",obj.multi=list(Mm=Mm,dat=Dd)))
})

 
          
       
setMethod("plot",signature(x="ALmult"), function(x,y=NULL,show.legend="right",...){
require(lattice)
load("GP.RData")                                                                                                                              #<<<<---- to be replaced by 'data(...)'
 dots <- list(...)
 sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
 if (is.null(dots$xlab)) dots$xlab <- "Length" ; if (is.null(dots$ylab)) dots$ylab <- "Age Proportion at Length" ; if (is.null(dots$main)) dots$main <- "Age-at-length multinomial analysis"
 
 
tab <- x@obj.multi$dat ; Mm <- x@obj.multi$Mm
strata <- colnames(tab)[ncol(tab)]
if (strata=="nostrata") show.legend <- ""
SMm <- summary(Mm,corr=FALSE)
Se <- SMm$standard.errors
Coef <- SMm$coefficients
ICm <- Coef-1.96*Se
ICp <- Coef+1.96*Se
tabAge <- tab[,1:(ncol(tab)-2)]
QijP <- tabAge/apply(tabAge,1,sum)
QijP[is.na(QijP)] <- 0
Q <- cbind.data.frame(as.data.frame(QijP),tab[,ncol(tab)-c(1:0)])
FIT <- Mm$fitted.values
test <- apply(tabAge,1,function(x) all(x==0))

QQ <- Q[!test,] ; L1 <- nrow(QQ) ; L2 <- ncol(QQ) ; N <- L2-2
datat <- data.frame(value=matrix(as.matrix(QQ[,1:N]),ncol=1),Age=rep(as.numeric(colnames(QQ)[1:N]),each=L1),Length=rep(QQ$Length,N),strata=rep(QQ[,strata],N))
datat$source <- "obs"
datatfit <- data.frame(value=matrix(as.matrix(FIT),ncol=1),Age=rep(as.numeric(colnames(QQ)[1:N]),each=L1),Length=rep(QQ$Length,N),strata=rep(QQ[,strata],N))
datatfit$source <- "smod"
datatfit2 <- rbind(datatfit,datat)      
  
ll <- nlevels(datatfit2$strata)      
xyplot(value ~ Length|Age, data = datatfit2, groups = interaction(strata,source), type = rep(c("p","l"),each=ll),
main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),
col=dots$l.col[1:ll],pch=dots$pch[1],cex=dots$p.cex[1],fill=dots$bg,lwd=dots$p.lwd[1],scales=list(font=dots$font.axis),distribute.type = TRUE,
key=list(points=list(pch=dots$pch[1],fill=dots$bg,cex=dots$p.cex[1],lwd=dots$p.lwd[1],col=dots$l.col[1:ll]),
                text=list(levels(datatfit2$strata)),title=strata,cex.title=0.8,space=show.legend,font=dots$font.lab,columns=1,border=TRUE))

})
       
 
 
########### 
# Example #
###########

load("sole2.RData")                                                                                                                  #<<<<---- to be replaced by 'data(...)'
 
aa <- AL.multi(sole2.cs,strata=NULL,age.plus=6)
bb1 <- AL.multi(sole2.cs,strata="quarter",age.plus=6)
bb2 <- AL.multi(sole2.cs,strata="quarter",elmts=c("2","3"),age.plus=8)      
       
                
plot(aa)       
plot(bb1,l.col=c("blue","red","gold"))       
plot(bb2,l.col=c("blue","red","gold"))


#Y a t'il moyen de comparer les différents échantillons d'une même strate?
#On a besoin de plus d'éléments en sortie pour conclure à des effets quelconque??
