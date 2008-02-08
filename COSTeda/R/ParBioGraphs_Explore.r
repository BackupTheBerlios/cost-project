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

data(GraphsPar)                                                                                                               
dots <- list(...)
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Length" ; if (is.null(dots$ylab)) dots$ylab <- "Weight" 
if (is.null(dots$main)) dots$main <- "Weight~Length" 

tab <- object@ca

if (selection){

print(xyplot(indWt~lenCls,data=tab,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),
      ylab=list(dots$ylab,font=dots$font.lab),pch=dots$pch[1],col=dots$col[1],cex=dots$p.cex[1],fill=dots$bg[1],lwd=dots$p.lwd[1],
      scales=list(font=dots$font.axis)))
trellis.focus("panel",1,1)
Reponse <- panel.identify()
return(list(l=Reponse,id.tab=tab[Reponse,]))

} else {

xyplot(indWt~lenCls,data=tab,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),
       ylab=list(dots$ylab,font=dots$font.lab),pch=dots$pch[1],col=dots$col[1],cex=dots$p.cex[1],fill=dots$p.bg[1],lwd=dots$p.lwd[1],
       scales=list(font=dots$font.axis))
}
})



#############################################
# scatterplot of individual maturity~length #
#############################################


setGeneric("ML.plot", function(object,...){
	standardGeneric("ML.plot")}
)



setMethod("ML.plot", signature(object="csData"), function(object,...){

data(GraphsPar)                                                                                                                  
dots <- list(...) 
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Length" ; if (is.null(dots$ylab)) dots$ylab <- "Maturity"
if (is.null(dots$main)) dots$main <- "Maturity~Length"
 
tab <- object@ca
tab$matStage <- factor(tab$matStage)

xyplot(matStage~lenCls,data=tab,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),
       ylab=list(dots$ylab,font=dots$font.lab),pch=dots$pch[1],col=dots$col[1],cex=dots$p.cex[1],fill=dots$p.bg[1],lwd=dots$p.lwd[1],
       scales=list(font=dots$font.axis))
})




#########################################
# boxplot of individual maturity~length #
#########################################



setGeneric("ML.boxplot", function(object,...){
	standardGeneric("ML.boxplot")}
)



setMethod("ML.boxplot", signature(object="csData"), function(object,...){

data(GraphsPar)                                                                                                           
dots <- list(...) ; if (is.null(dots$pch)) dots$pch <- "||"
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Length" ; if (is.null(dots$ylab)) dots$ylab <- "Maturity"
if (is.null(dots$main)) dots$main <- "Maturity~Length"

tab <- object@ca
tab$matStage <- factor(tab$matStage)

bwplot(matStage~lenCls,data=tab,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),
       ylab=list(dots$ylab,font=dots$font.lab),pch=dots$pch[1],cex=2,fill=dots$p.bg[1],scales=list(font=dots$font.axis),
       par.settings=list(box.rectangle=list(col=dots$col[1]),box.umbrella=list(col=dots$col[1],lty=dots$lty[1]),
       plot.symbol=list(col=dots$col[1])))
})




####################################
# boxplot of individual sex~length #
####################################



setGeneric("SL.boxplot", function(object,...){
	standardGeneric("SL.boxplot")}
)


setMethod("SL.boxplot", signature(object="csData"), function(object,main="",xlab="Length",ylab="Sex",...){

data(GraphsPar)                                                                                                                  
dots <- list(...) ; if (is.null(dots$pch)) dots$pch <- "||"
sapply(names(GP),function(x) if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) dots$xlab <- "Length" ; if (is.null(dots$ylab)) dots$ylab <- "Sex"
if (is.null(dots$main)) dots$main <- "Sex~Length"

tab <- object@ca
tab$sex <- factor(as.character(tab$sex),exclude="U")

bwplot(sex~lenCls,data=tab,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),
       ylab=list(dots$ylab,font=dots$font.lab),pch=dots$pch[1],cex=2,fill=dots$p.bg[1],scales=list(font=dots$font.axis),
       par.settings=list(box.rectangle=list(col=dots$col[1]),box.umbrella=list(col=dots$col[1],lty=dots$lty[1]),
       plot.symbol=list(col=dots$col[1])))
})





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

