###########################################
##                                       ##
## Plots of Other biological  parameters ##
##                                       ##
##            MM 07/02/2008              ##
###########################################


###########################################
# scatterplot of individual weight~length #
###########################################

setGeneric("wlPlot", function(object,
                              selection=FALSE,...){
	standardGeneric("wlPlot")}
)




setMethod("wlPlot", signature(object="csData"), function(object,
                                                         selection=FALSE,...){

tab <- object@ca 

#tests on ca table fields
if (all(is.na(tab$indWt))) stop("no individual weight data in ca table!!")
if (all(is.na(tab$lenCls))) stop("no length class data in ca table!!") 

data(GraphsPar,envir=environment())                                                                                                               
dots <- list(...)
sapply(names(GP),function(x) 
                  if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) 
                    eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) 
  dots$xlab <- "Length(mm)"
if (is.null(dots$ylab)) 
  dots$ylab <- "Weight(g)" 
if (is.null(dots$main)) 
  dots$main <- "Scatter plot of individual weight at length" 

#missing length classes are taken into account
lenC <- c(1,5,10,25)
names(lenC) <- c("mm","scm","cm","25mm")
tab$lenCls <- factor(tab$lenCls,levels=seq(min(tab$lenCls),max(tab$lenCls),by=lenC[as.character(tab$lenCode[1])]))

if (selection){

  print(xyplot(indWt~lenCls,data=tab,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),
               ylab=list(dots$ylab,font=dots$font.lab),pch=dots$pch[1],col=dots$col[1],cex=dots$p.cex[1],fill=dots$p.bg[1],lwd=dots$p.lwd[1],
               scales=list(font=dots$font.axis,x=list(rot=dots$rot[1])),drop.unused.levels=FALSE))
  trellis.focus("panel",1,1)
  Reponse <- panel.identify()
  id.tab <- tab[Reponse,] 
  tabOcc <- paste(tab$indWt,tab$lenCls,sep=":::")
  idOcc <- paste(id.tab$indWt,id.tab$lenCls,sep=":::")
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


setGeneric("mlPlot", function(object,
                              selection=FALSE,
                              ...){
	standardGeneric("mlPlot")}
)



setMethod("mlPlot", signature(object="csData"), function(object,
                                                         selection=FALSE,
                                                         ...){
tab <- object@ca 

#tests on ca table fields
if (all(is.na(tab$matStage))) stop("no maturity stage data in ca table!!")
if (all(is.na(tab$lenCls))) stop("no length class data in ca table!!") 


data(GraphsPar,envir=environment())                                                                                                                  
dots <- list(...) 
sapply(names(GP),function(x) 
                  if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) 
                    eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) 
  dots$xlab <- "Length(mm)" 
if (is.null(dots$ylab)) 
  dots$ylab <- "Maturity"
if (is.null(dots$main)) 
  dots$main <- "Scatter plot of individual maturity at length"
 
tab$matStage <- factor(tab$matStage)
#missing length classes are taken into account
lenC <- c(1,5,10,25)
names(lenC) <- c("mm","scm","cm","25mm")
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


setGeneric("slPlot", function(object,
                              selection=FALSE,
                              ...){
	standardGeneric("slPlot")}
)



setMethod("slPlot", signature(object="csData"), function(object,
                                                         selection=FALSE,
                                                         ...){

tab <- object@ca 

#tests on ca table fields
if (all(is.na(tab$sex))) stop("no sex data in ca table!!")
if (all(is.na(tab$lenCls))) stop("no length class data in ca table!!") 


data(GraphsPar,envir=environment())                                                                                                                  
dots <- list(...) 
sapply(names(GP),function(x) 
                  if (is.null(eval(parse('',text=paste("dots$",x,sep="")))))
                    eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) 
  dots$xlab <- "Length(mm)"
if (is.null(dots$ylab)) 
  dots$ylab <- "Sex"
if (is.null(dots$main))
  dots$main <- "Scatter plot of individual sex at length"
 
tab$sex <- factor(as.character(tab$sex),exclude="U")
#missing length classes are taken into account
lenC <- c(1,5,10,25)
names(lenC) <- c("mm","scm","cm","25mm")
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



setGeneric("wlBoxplot", function(object,...){
	standardGeneric("wlBoxplot")}
)



setMethod("wlBoxplot", signature(object="csData"), function(object,...){

tab <- object@ca 

#tests on ca table fields
if (all(is.na(tab$indWt))) stop("no individual weight data in ca table!!")
if (all(is.na(tab$lenCls))) stop("no length class data in ca table!!") 

data(GraphsPar,envir=environment())                                                                                                           
dots <- list(...)
if (is.null(dots$pch)) dots$pch <- 20
sapply(names(GP),function(x) 
                  if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) 
                    eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) 
  dots$xlab <- "Length(mm)" 
if (is.null(dots$ylab)) 
  dots$ylab <- "Weight(g)"
if (is.null(dots$main)) 
  dots$main <- "Boxplot of individual weight at length"

#missing length classes are taken into account
lenC <- c(1,5,10,25)
names(lenC) <- c("mm","scm","cm","25mm")
tab$lenCls <- factor(tab$lenCls,levels=seq(min(tab$lenCls),max(tab$lenCls),by=lenC[as.character(tab$lenCode[1])]))

bwplot(indWt~lenCls,data=tab,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),
       ylab=list(dots$ylab,font=dots$font.lab),pch=dots$pch[1],cex=1.6,fill=dots$p.bg[1],scales=list(font=dots$font.axis,x=list(rot=dots$rot[1])),
       par.settings=list(box.rectangle=list(col=dots$col[1]),box.umbrella=list(col=dots$col[1],lty=dots$lty[1]),
       plot.symbol=list(col=dots$col[1])),drop.unused.levels=FALSE)
})





#########################################
# boxplot of individual maturity~length #
#########################################



setGeneric("mlBoxplot", function(object,...){
	standardGeneric("mlBoxplot")}
)



setMethod("mlBoxplot", signature(object="csData"), function(object,...){

tab <- object@ca

#tests on ca table fields
if (all(is.na(tab$matStage))) stop("no maturity stage data in ca table!!")
if (all(is.na(tab$lenCls))) stop("no length class data in ca table!!") 


data(GraphsPar,envir=environment())                                                                                                           
dots <- list(...)
if (is.null(dots$pch)) dots$pch <- 20
sapply(names(GP),function(x) 
                  if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) 
                    eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab)) 
  dots$xlab <- "Length(mm)"
if (is.null(dots$ylab)) 
  dots$ylab <- "Maturity"
if (is.null(dots$main)) 
  dots$main <- "Boxplot of individual maturity at length"

tab$matStage <- factor(tab$matStage)
#missing length classes are taken into account
lenC <- c(1,5,10,25)
names(lenC) <- c("mm","scm","cm","25mm")
tab$lenCls <- factor(tab$lenCls,levels=seq(min(tab$lenCls),max(tab$lenCls),by=lenC[as.character(tab$lenCode[1])]))

bwplot(matStage~lenCls,data=tab,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),
       ylab=list(dots$ylab,font=dots$font.lab),pch=dots$pch[1],cex=2,fill=dots$p.bg[1],scales=list(font=dots$font.axis,x=list(rot=dots$rot[1])),
       par.settings=list(box.rectangle=list(col=dots$col[1]),box.umbrella=list(col=dots$col[1],lty=dots$lty[1]),
       plot.symbol=list(col=dots$col[1])),drop.unused.levels=FALSE)
})




####################################
# boxplot of individual sex~length #
####################################



setGeneric("slBoxplot", function(object,...){
	standardGeneric("slBoxplot")}
)


setMethod("slBoxplot", signature(object="csData"), function(object,...){

tab <- object@ca

#tests on ca table fields
if (all(is.na(tab$sex))) stop("no sex data in ca table!!")
if (all(is.na(tab$lenCls))) stop("no length class data in ca table!!") 

data(GraphsPar,envir=environment())                                                                                                                  
dots <- list(...)
if (is.null(dots$pch)) dots$pch <- 20
sapply(names(GP),function(x) 
                  if (is.null(eval(parse('',text=paste("dots$",x,sep=""))))) 
                    eval(parse('',text=paste("dots$",x," <<- GP$",x,sep=""))))
if (is.null(dots$xlab))
  dots$xlab <- "Length(mm)"
if (is.null(dots$ylab)) 
  dots$ylab <- "Sex"
if (is.null(dots$main)) 
  dots$main <- "Boxplot of individual sex at length"

tab$sex <- factor(as.character(tab$sex),exclude="U")
#missing length classes are taken into account
lenC <- c(1,5,10,25)
names(lenC) <- c("mm","scm","cm","25mm")
tab$lenCls <- factor(tab$lenCls,levels=seq(min(tab$lenCls),max(tab$lenCls),by=lenC[as.character(tab$lenCode[1])]))

bwplot(sex~lenCls,data=tab,main=list(dots$main,font=dots$font.main),xlab=list(dots$xlab,font=dots$font.lab),
       ylab=list(dots$ylab,font=dots$font.lab),pch=dots$pch[1],cex=2,fill=dots$p.bg[1],scales=list(font=dots$font.axis,x=list(rot=dots$rot[1])),
       par.settings=list(box.rectangle=list(col=dots$col[1]),box.umbrella=list(col=dots$col[1],lty=dots$lty[1]),
       plot.symbol=list(col=dots$col[1])),drop.unused.levels=FALSE)
})


##############################################################################
##############################################################################
######################### Methods to export ##################################
##############################################################################
##############################################################################


########################
# Specific plot.design #
########################

setGeneric("csPlot.design", function(object,...){
	standardGeneric("csPlot.design")}
)


setMethod("csPlot.design", signature(object="csDataVal"), function(object,...){
 
tab <- ca(object)
tab$age <- factor(tab$age,exclude=NA)
tab$matStage <- factor(tab$matStage,exclude=NA)
tab$sex <- factor(as.character(tab$sex),exclude="U")
tab$lenCls <- as.numeric(as.character(tab$lenCls))

des.tab <- tab[,c("age","matStage","sex","lenCls","indWt")]
#index of 'empty' columns
test <- apply(des.tab,2,function(x) all(is.na(x)))
des.tab <- des.tab[,!test]

if (ncol(des.tab)==0) stop("no data to be used in object!!")

plot.design(des.tab,...) 
})



#################
# Bio Par plots #
#################


setGeneric("bioPar.plot", function(object,
                                   type="wl",     #or "ml" or "sl"
                                   selection=FALSE,...){
	standardGeneric("bioPar.plot")}
)


setMethod("bioPar.plot", signature(object="csData"), function(object,
                                                              type="wl",
                                                              selection=FALSE,
                                                              ...){
eval(parse('',text=paste(type,"Plot(object,selection=selection,...)",sep="")))
})


setGeneric("bioPar.boxplot", function(object,
                                      type="wl",     #or "ml" or "sl"
                                      ...){
	standardGeneric("bioPar.boxplot")}
)


setMethod("bioPar.boxplot", signature(object="csData"), function(object,
                                                                 type="wl",
                                                                 ...){
eval(parse('',text=paste(type,"Boxplot(object,...)",sep="")))
})




