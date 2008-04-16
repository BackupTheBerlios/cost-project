################################################
##                                            ##
## Age at Length Analysis (multinomial model) ##
##                                            ##
##       MM 04/02/2008                        ##
################################################


setClass("ageLenMult",representation(objMult="list"),
                      prototype(objMult=list()))

 
setGeneric("ageLenMultinom", function(object,
                                timeStrata=NULL, #ex: "year","quarter","month"
                                spaceStrata=NULL, #ex: "area","rect"
                                techStrata=NULL, #ex: "commCat"
                                elmts=list(tp="all",sp="all",tc="all"),
                                grps=NULL, #ex: "timeStrata","spaceStrata","techStrata"
                                age.plus=-1,
                                  ...){
	standardGeneric("ageLenMultinom")}
)


#no options for multinomial modelisation for the moment  ----> to be continued
       

setMethod("ageLenMultinom", signature(object="csDataVal"), function(object,
                                                                timeStrata=NULL, #ex: "year","quarter","month"
                                                                spaceStrata=NULL, #ex: "area","rect"
                                                                techStrata=NULL, #ex: "commCat"
                                                                elmts=list(tp="all",sp="all",tc="all"),
                                                                grps=NULL, #ex: "timeStrata","spaceStrata","techStrata"
                                                                age.plus=-1, #ex: 6, 8
                                                                ...){

require(nnet)
tab <- object@ca
tab$age <- as.numeric(as.character(tab$age)) 

#test on elmts
if (!all(names(elmts)%in%c("tp","sp","tc"))) stop("'elmts' parameter is not well defined") 
#grouping ages
if (age.plus>=0) tab$age[tab$age>=age.plus] <- age.plus
newtab <- tab[,c("age","lenCls",timeStrata,spaceStrata,techStrata)] 

#subsetting
tpT <- !is.null(timeStrata)
if (tpT&(!("all"%in%elmts$tp))) newtab <- newtab[as.character(newtab[,timeStrata])%in%elmts$tp,]

spT <- !is.null(spaceStrata) 
if (spT&(!("all"%in%elmts$sp))) newtab <- newtab[as.character(newtab[,spaceStrata])%in%elmts$sp,]

tcT <- !is.null(techStrata) 
if (tcT&(!("all"%in%elmts$tc))) newtab <- newtab[as.character(newtab[,techStrata])%in%elmts$tc,]


if (tpT) newtab[,timeStrata] <- factor(newtab[,timeStrata]) 
if (spT) newtab[,spaceStrata] <- factor(newtab[,spaceStrata]) 
if (tcT) newtab[,techStrata] <- factor(newtab[,techStrata])

alk <- tapply(newtab$age,lapply(c("lenCls","age",timeStrata,spaceStrata,techStrata),function(x) newtab[,x]),length)
alk[is.na(alk)] <- 0 
alkNam <- dimnames(alk)
alkNam[[2]] <- NULL
matTest <- apply(alk,seq(from=1,to=length(c(timeStrata,spaceStrata,techStrata))+2)[-2],function(x) all(x==0))
test <- matrix(matTest,ncol=1)
Age <- apply(alk,2,rbind)
Dim <- dim(alk)[-2]
Ld <- length(Dim)

Tab2 <- do.call("cbind",lapply(1:Ld,function(x) {if ((x+1)>Ld) 
                                                  vec1 <- 1 
                                                 else 
                                                  vec1 <- Dim[(x+1):Ld]
                                                 if (x==1) 
                                                  vec2 <- 1 
                                                 else 
                                                  vec2 <- Dim[1:(x-1)]
                                                 return(rep(rep(alkNam[[x]],prod(vec1)),each=prod(vec2)))}))

Tab2 <- as.data.frame(Tab2)
names(Tab2) <- c("Length",timeStrata,spaceStrata,techStrata)
Tab2$Length <- as.numeric(as.character(Tab2$Length))
Tab <- cbind.data.frame(as.data.frame(Age),Tab2)
                          
TAB <- Tab[!test,] 
Age <- Age[!test,]
ntest <- all(is.null(timeStrata),is.null(spaceStrata),is.null(techStrata))
eval(parse('',text=paste("Mm <- multinom(Age~Length",c(paste("*",paste(c(timeStrata,spaceStrata,techStrata),collapse="*"),"-1",sep=""),"")[c(!ntest,ntest)],
                         ",Hess=T,model=T,data=TAB)",sep="")))
                   
return(new("ageLenMult",objMult=list(timeStrata=timeStrata,
                                     spaceStrata=spaceStrata,
                                     techStrata=techStrata,
                                     grps=grps,
                                     Mm=Mm,
                                     dat=Tab)))
})

 
          
       
setMethod("plot",signature(x="ageLenMult"), function(x, #'ageLenMult' object
                                                     y=NULL,
                                                     show.legend="right", #ex: "", "left"
                                                     ...){

data(GraphsPar)                                                                                                                            
dots <- list(...)
sapply(names(GP),function(x) 
                  if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) 
                    eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) 
  dots$xlab <- "Length (mm)"
if (is.null(dots$ylab)) 
  dots$ylab <- "Age Proportion at Length" 
if (is.null(dots$main)) 
  dots$main <- "Age-at-length multinomial analysis"
 
timeStrata <- x@objMult$timeStrata
spaceStrata <- x@objMult$spaceStrata 
techStrata <- x@objMult$techStrata  
tab <- x@objMult$dat
Mm <- x@objMult$Mm
stratas <- c(x@objMult$timeStrata,x@objMult$spaceStrata,x@objMult$techStrata)

Ls <- length(stratas)
if (Ls==0) show.legend <- ""
SMm <- summary(Mm,corr=FALSE)
lim <- grep("Length",names(tab))
tabAge <- tab[,1:(lim-1)]
QijP <- tabAge/apply(tabAge,1,sum)
QijP[is.na(QijP)] <- 0
Q <- cbind.data.frame(as.data.frame(QijP),tab[,lim:ncol(tab),drop=FALSE])
FIT <- Mm$fitted.values
test <- apply(tabAge,1,function(x) all(x==0))

QQ <- Q[!test,]
L1 <- nrow(QQ)
L2 <- ncol(QQ) 
eval(parse('',text=paste("datat <- data.frame(value=matrix(as.matrix(QQ[,1:(lim-1)]),ncol=1),Age=rep(as.numeric(colnames(QQ)[1:(lim-1)]),each=L1),Length=rep(QQ$Length,(lim-1))",
                         ","[Ls!=0],paste(paste(stratas,"=rep(QQ$",stratas,",(lim-1))",sep=""),collapse=",")[Ls!=0],")",sep="")))
datat$source <- "obs"
eval(parse('',text=paste("datatfit <- data.frame(value=matrix(as.matrix(FIT),ncol=1),Age=rep(as.numeric(colnames(QQ)[1:(lim-1)]),each=L1),Length=rep(QQ$Length,(lim-1))",
                         ","[Ls!=0],paste(paste(stratas,"=rep(QQ$",stratas,",(lim-1))",sep=""),collapse=",")[Ls!=0],")",sep="")))
datatfit$source <- "smod"
datatfit2 <- rbind(datatfit,datat)      
datatfit2$Age <- factor(as.character(datatfit2$Age),levels=as.character(sort(unique(as.numeric(as.character(datatfit2$Age))))))

if (is.null(x@objMult$grps)) 
  grps <- NULL 
else 
  grps <- eval(parse('',text=x@objMult$grps))  
if (Ls==0) grps <- NULL

if (is.null(grps)) {
  grps <- "nostrata" 
  datatfit2$nostrata <- factor(rep("all",nrow(datatfit2)))}

eval(parse('',text=paste("ll <- nlevels(datatfit2$",grps,")",sep="")))      

declStr <- stratas[!stratas%in%grps]
#erase the last page if length(declStr)>=2
lll <- length(declStr)
if (lll>1) 
  pageNb <- sum(sapply(declStr[2:lll],function(x) eval(parse('',text=paste("length(unique(QQ$",x,"))",sep=""))))) 
else 
  pageNb <- 1 

strip.col <- trellis.par.get("strip.background")$col
gpTest <- grps=="nostrata"
 
eval(parse('',text=paste("xyplot(value~Length|",paste(c("Age",declStr),collapse="*"),",data = datatfit2, groups = interaction(",grps,",source), type = rep(c(\"p\",\"l\"),each=ll),",
"main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),ylab=list(dots$ylab,font=dots$font.lab),par.strip.text=list(font=dots$font.lab),",
"col=dots$l.col[1:ll],pch=dots$pch[1],cex=dots$p.cex[1],fill=dots$bg,lwd=dots$p.lwd[1],scales=list(font=dots$font.axis),distribute.type = TRUE,page=function(n) {dev.copy(x11) ; if (n==pageNb) dev.off()},",
"key=list(points=list(pch=c(rep(dots$pch[1],ll),NA,rep(15,lll+1)),fill=dots$bg,cex=dots$p.cex[1],lwd=dots$p.lwd[1],col=c(dots$l.col[1:ll],NA,strip.col[1:(lll+1)])),",
"text=list(c(levels(datatfit2$",grps,"),\"\",\"age\",declStr)),title=",c(paste("\"",grps,"\"",sep=""),"NULL")[c(!gpTest,gpTest)],
",cex.title=0.8,space=show.legend,font=dots$font.lab,columns=1,border=TRUE))",sep="")))

})
  
       
       
setGeneric("plot.ageLenMult", function(x, 
                                       show.legend="right", 
                                       ...){
	standardGeneric("plot.ageLenMult")}
)



setMethod("plot.ageLenMult",signature(x="ageLenMult"), function(x, 
                                                                show.legend="right",
                                                                ...){
plot(x=x,show.legend=show.legend,...)
})


