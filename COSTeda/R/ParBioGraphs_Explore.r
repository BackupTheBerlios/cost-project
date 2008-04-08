###########################################
##                                       ##
## Plots of Other biological  parameters ##
##                                       ##
##            MM 07/02/2008              ##
###########################################


###########################################
# scatterplot of individual weight~length #
###########################################

setGeneric("WL.plot", function(object,
                               selection=FALSE,...){
	standardGeneric("WL.plot")}
)




setMethod("WL.plot", signature(object="csData"), function(object,
                                                          selection=FALSE,...){
tab <- object@ca 

data(GraphsPar)                                                                                                               
dots <- list(...)
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Length(mm)" ; if (is.null(dots$ylab)) dots$ylab <- "Weight(g)" 
if (is.null(dots$main)) dots$main <- "Scatter plot of individual weight at length" 

#on recode les classes de tailles pour combler les classes manquantes
lenC <- c(1,5,10,25) ; names(lenC) <- c("mm","scm","cm","25mm")
tab$lenCls <- factor(tab$lenCls,levels=seq(min(tab$lenCls),max(tab$lenCls),by=lenC[as.character(tab$lenCode[1])]))

if (selection){

print(xyplot(indWt~lenCls,data=tab,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),
      ylab=list(dots$ylab,font=dots$font.lab),pch=dots$pch[1],col=dots$col[1],cex=dots$p.cex[1],fill=dots$p.bg[1],lwd=dots$p.lwd[1],
      scales=list(font=dots$font.axis,x=list(rot=dots$rot[1])),drop.unused.levels=FALSE))
trellis.focus("panel",1,1)
Reponse <- panel.identify()
id.tab <- tab[Reponse,] 
tabOcc <- paste(tab$indWt,tab$lenCls,sep=":::") ; idOcc <- paste(id.tab$indWt,id.tab$lenCls,sep=":::")
invisible(list(l=Reponse,id.tab=tab[tabOcc%in%idOcc,]))

} else {

xyplot(indWt~lenCls,data=tab,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),
       ylab=list(dots$ylab,font=dots$font.lab),pch=dots$pch[1],col=dots$col[1],cex=dots$p.cex[1],fill=dots$p.bg[1],lwd=dots$p.lwd[1],
       scales=list(font=dots$font.axis,x=list(rot=dots$rot[1])),drop.unused.levels=FALSE)
}
})



#############################################
# scatterplot of individual maturity~length #
#############################################


setGeneric("ML.plot", function(object,
                               selection=FALSE,
                               ...){
	standardGeneric("ML.plot")}
)



setMethod("ML.plot", signature(object="csData"), function(object,
                                                          selection=FALSE,
                                                          ...){
tab <- object@ca 

data(GraphsPar)                                                                                                                  
dots <- list(...) 
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Length(mm)" ; if (is.null(dots$ylab)) dots$ylab <- "Maturity"
if (is.null(dots$main)) dots$main <- "Scatter plot of individual maturity at length"
 
tab$matStage <- factor(tab$matStage)
#on recode les classes de tailles pour combler les classes manquantes
lenC <- c(1,5,10,25) ; names(lenC) <- c("mm","scm","cm","25mm")
tab$lenCls <- factor(tab$lenCls,levels=seq(min(tab$lenCls),max(tab$lenCls),by=lenC[as.character(tab$lenCode[1])]))


if (selection){

print(xyplot(matStage~lenCls,data=tab,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),
      ylab=list(dots$ylab,font=dots$font.lab),pch=dots$pch[1],col=dots$col[1],cex=dots$p.cex[1],fill=dots$p.bg[1],lwd=dots$p.lwd[1],
      scales=list(font=dots$font.axis,x=list(rot=dots$rot[1])),drop.unused.levels=FALSE))
trellis.focus("panel",1,1)
Reponse <- panel.identify()
id.tab <- tab[Reponse,] 
tabOcc <- paste(tab$matStage,tab$lenCls,sep=":::") ; idOcc <- paste(id.tab$matStage,id.tab$lenCls,sep=":::")
invisible(list(l=Reponse,id.tab=tab[tabOcc%in%idOcc,]))

} else {

xyplot(matStage~lenCls,data=tab,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),
       ylab=list(dots$ylab,font=dots$font.lab),pch=dots$pch[1],col=dots$col[1],cex=dots$p.cex[1],fill=dots$p.bg[1],lwd=dots$p.lwd[1],
       scales=list(font=dots$font.axis,x=list(rot=dots$rot[1])),drop.unused.levels=FALSE)
}
})



########################################
# scatterplot of individual sex~length #
########################################


setGeneric("SL.plot", function(object,
                               selection=FALSE,
                               ...){
	standardGeneric("SL.plot")}
)



setMethod("SL.plot", signature(object="csData"), function(object,
                                                          selection=FALSE,
                                                          ...){

tab <- object@ca 

data(GraphsPar)                                                                                                                  
dots <- list(...) 
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Length(mm)" ; if (is.null(dots$ylab)) dots$ylab <- "Sex"
if (is.null(dots$main)) dots$main <- "Scatter plot of individual sex at length"
 
tab$sex <- factor(as.character(tab$sex),exclude="U")
#on recode les classes de tailles pour combler les classes manquantes
lenC <- c(1,5,10,25) ; names(lenC) <- c("mm","scm","cm","25mm")
tab$lenCls <- factor(tab$lenCls,levels=seq(min(tab$lenCls),max(tab$lenCls),by=lenC[as.character(tab$lenCode[1])]))

if (selection){

print(xyplot(sex~lenCls,data=tab,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),
      ylab=list(dots$ylab,font=dots$font.lab),pch=dots$pch[1],col=dots$col[1],cex=dots$p.cex[1],fill=dots$p.bg[1],lwd=dots$p.lwd[1],
      scales=list(font=dots$font.axis,x=list(rot=dots$rot[1])),drop.unused.levels=FALSE))
trellis.focus("panel",1,1)
Reponse <- panel.identify()
id.tab <- tab[Reponse,] 
tabOcc <- paste(tab$sex,tab$lenCls,sep=":::") ; idOcc <- paste(id.tab$sex,id.tab$lenCls,sep=":::")
invisible(list(l=Reponse,id.tab=tab[tabOcc%in%idOcc,]))

} else {

xyplot(sex~lenCls,data=tab,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),
       ylab=list(dots$ylab,font=dots$font.lab),pch=dots$pch[1],col=dots$col[1],cex=dots$p.cex[1],fill=dots$p.bg[1],lwd=dots$p.lwd[1],
       scales=list(font=dots$font.axis,x=list(rot=dots$rot[1])),drop.unused.levels=FALSE)
}
})



#######################################
# boxplot of individual weight~length #
#######################################



setGeneric("WL.boxplot", function(object,...){
	standardGeneric("WL.boxplot")}
)



setMethod("WL.boxplot", signature(object="csData"), function(object,...){

tab <- object@ca 

data(GraphsPar)                                                                                                           
dots <- list(...) ; if (is.null(dots$pch)) dots$pch <- 20
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Length(mm)" ; if (is.null(dots$ylab)) dots$ylab <- "Weight(g)"
if (is.null(dots$main)) dots$main <- "Boxplot of individual weight at length"

#on recode les classes de tailles pour combler les classes manquantes
lenC <- c(1,5,10,25) ; names(lenC) <- c("mm","scm","cm","25mm")
tab$lenCls <- factor(tab$lenCls,levels=seq(min(tab$lenCls),max(tab$lenCls),by=lenC[as.character(tab$lenCode[1])]))

bwplot(indWt~lenCls,data=tab,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),
       ylab=list(dots$ylab,font=dots$font.lab),pch=dots$pch[1],cex=1.6,fill=dots$p.bg[1],scales=list(font=dots$font.axis,x=list(rot=dots$rot[1])),
       par.settings=list(box.rectangle=list(col=dots$col[1]),box.umbrella=list(col=dots$col[1],lty=dots$lty[1]),
       plot.symbol=list(col=dots$col[1])),drop.unused.levels=FALSE)
})





#########################################
# boxplot of individual maturity~length #
#########################################



setGeneric("ML.boxplot", function(object,...){
	standardGeneric("ML.boxplot")}
)



setMethod("ML.boxplot", signature(object="csData"), function(object,...){

tab <- object@ca

data(GraphsPar)                                                                                                           
dots <- list(...) ; if (is.null(dots$pch)) dots$pch <- 20
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Length(mm)" ; if (is.null(dots$ylab)) dots$ylab <- "Maturity"
if (is.null(dots$main)) dots$main <- "Boxplot of individual maturity at length"

tab$matStage <- factor(tab$matStage)
#on recode les classes de tailles pour combler les classes manquantes
lenC <- c(1,5,10,25) ; names(lenC) <- c("mm","scm","cm","25mm")
tab$lenCls <- factor(tab$lenCls,levels=seq(min(tab$lenCls),max(tab$lenCls),by=lenC[as.character(tab$lenCode[1])]))

bwplot(matStage~lenCls,data=tab,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),
       ylab=list(dots$ylab,font=dots$font.lab),pch=dots$pch[1],cex=2,fill=dots$p.bg[1],scales=list(font=dots$font.axis,x=list(rot=dots$rot[1])),
       par.settings=list(box.rectangle=list(col=dots$col[1]),box.umbrella=list(col=dots$col[1],lty=dots$lty[1]),
       plot.symbol=list(col=dots$col[1])),drop.unused.levels=FALSE)
})




####################################
# boxplot of individual sex~length #
####################################



setGeneric("SL.boxplot", function(object,...){
	standardGeneric("SL.boxplot")}
)


setMethod("SL.boxplot", signature(object="csData"), function(object,main="",xlab="Length",ylab="Sex",...){

tab <- object@ca

data(GraphsPar)                                                                                                                  
dots <- list(...) ; if (is.null(dots$pch)) dots$pch <- 20
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Length(mm)" ; if (is.null(dots$ylab)) dots$ylab <- "Sex"
if (is.null(dots$main)) dots$main <- "Boxplot of individual sex at length"

tab$sex <- factor(as.character(tab$sex),exclude="U")
#on recode les classes de tailles pour combler les classes manquantes
lenC <- c(1,5,10,25) ; names(lenC) <- c("mm","scm","cm","25mm")
tab$lenCls <- factor(tab$lenCls,levels=seq(min(tab$lenCls),max(tab$lenCls),by=lenC[as.character(tab$lenCode[1])]))

bwplot(sex~lenCls,data=tab,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),
       ylab=list(dots$ylab,font=dots$font.lab),pch=dots$pch[1],cex=2,fill=dots$p.bg[1],scales=list(font=dots$font.axis,x=list(rot=dots$rot[1])),
       par.settings=list(box.rectangle=list(col=dots$col[1]),box.umbrella=list(col=dots$col[1],lty=dots$lty[1]),
       plot.symbol=list(col=dots$col[1])),drop.unused.levels=FALSE)
})


################################################################################
################################################################################
######################### Méthodes à exporter ##################################
################################################################################
################################################################################


########################
# Specific plot.design #
########################

setGeneric("L.plot.design", function(object,...){
	standardGeneric("L.plot.design")}
)


setMethod("L.plot.design", signature(object="csData"), function(object,...){
 
tab <- object@ca
tab$age <- factor(tab$age,exclude=NA) ; tab$matStage <- factor(tab$matStage,exclude=NA) ; tab$sex <- factor(as.character(tab$sex),exclude="U")
des.tab <- tab[,c("age","matStage","sex","lenCls","indWt")]

plot.design(des.tab,...) 
})



#################
# Bio Par plots #
#################


setGeneric("BioPar.plot", function(object,
                                   type="WL",     #ou "ML" ou "SL"
                                   selection=FALSE,...){
	standardGeneric("BioPar.plot")}
)


setMethod("BioPar.plot", signature(object="csData"), function(object,
                                                              type="WL",
                                                              selection=FALSE,
                                                              ...){
eval(parse('',text=paste(type,".plot(object,selection=selection,...)",sep="")))
})


setGeneric("BioPar.boxplot", function(object,
                                      type="WL",     #ou "ML" ou "SL"
                                      ...){
	standardGeneric("BioPar.boxplot")}
)


setMethod("BioPar.boxplot", signature(object="csData"), function(object,
                                                                 type="WL",
                                                                 ...){
eval(parse('',text=paste(type,".boxplot(object,...)",sep="")))
})




