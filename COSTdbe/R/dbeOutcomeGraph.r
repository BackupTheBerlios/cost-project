

#library(COSTcore)
#library(lattice)
#source("C:/Documents and Settings/mmerzere/Bureau/outcomeObject.r")
# 

#############################################
#                                           #
# Graphical display from 'dbeOutput' object #
#                                           #
#############################################

        
setGeneric("dbePlot", function(object,
                               Slot,                                 
                               type="bar",                           
                               Xstratum=NULL,
                               step=NA,                        
                               dispKey=TRUE,                         
                               ...){
standardGeneric("dbePlot")}
)



 
setMethod("dbePlot",signature(object="dbeOutput"),        
           function(object,                               #'dbeOutput' object 
                    Slot,                                 #ex: "lenStruc"
                    type="bar",                           #to be chosen between "bar", "line", "point" 
                    Xstratum=NULL,                        #stratum displayed on x-axis if Slot is in c("nSamp","nMes","totalN","totalNvar","totalW","totalWvar")
                    step=NA,                              #length or age classes step (if NA, empty classes are not drawn)
                    dispKey=TRUE,                         #if TRUE and if various panels are displayed, a describing key is displayed
                    ...){                                 #further graphical parameters

if (!Slot%in%slotNames(object)[7:16]) stop("Wrong 'Slot' parameter!!")
if (!type%in%c("bar","point","line")) stop("Wrong 'type' parameter!!")

  #-----------------------------------------------------------------------------
  # Graphical parameters
  #-----------------------------------------------------------------------------

dots <- list(...)  

if (is.null(dots$col)) dots$col <- "skyblue"
if (is.null(dots$lwd)) dots$lwd <- 1
if (is.null(dots$lty)) dots$lty <- 1
if (is.null(dots$pch)) dots$pch <- 21

if (is.null(dots$cex)) dots$cex <- 0.8
if (is.null(dots$p.cex)) dots$p.cex <- 1.2
if (is.null(dots$k.cex)) dots$k.cex <- 1.2
if (is.null(dots$cex.lab)) dots$cex.lab <- 1.1
if (is.null(dots$cex.sub)) dots$cex.sub <- 1.1
if (is.null(dots$cex.axis)) dots$cex.axis <- 1
if (is.null(dots$cex.main)) dots$cex.main <- 1.2

if (is.null(dots$col.lab)) dots$col.lab <- "black"
if (is.null(dots$col.sub)) dots$col.sub <- "black"
if (is.null(dots$col.axis)) dots$col.axis <- "black"
if (is.null(dots$col.main)) dots$col.main <- "black"

if (is.null(dots$font)) dots$font <- 6
if (is.null(dots$font.lab)) dots$font.lab <- 7
if (is.null(dots$font.sub)) dots$font.sub <- 6
if (is.null(dots$font.axis)) dots$font.axis <- 7
if (is.null(dots$font.main)) dots$font.main <- 7

if (is.null(dots$rot)) dots$rot <- 90

  #-----------------------------------------------------------------------------
  # Extraction of the numerical data, and formatting process
  #-----------------------------------------------------------------------------

tab <- slot(object,Slot) ; if (length(tab)<3) tab <- tab$estim  #if slot is a list with 2 elements
if (nrow(tab)==0) stop("no available data!!")
if (all(is.na(tab))) stop("no available data!!")
if (all(levels(factor(as.character(tab$time)))=="all")) tab$time <- NA
if (all(levels(factor(as.character(tab$space)))=="all")) tab$space <- NA
if (all(levels(factor(as.character(tab$technical)))=="all")) tab$technical <- NA

timeStrata <- spaceStrata <- techStrata <- TRUE
if (all(is.na(tab$time))) timeStrata <- FALSE
if (all(is.na(tab$space))) spaceStrata <- FALSE
if (all(is.na(tab$technical))) techStrata <- FALSE

#indicates if Slot is in c("nSamp","nMes","totalN","totalNvar","totalW","totalWvar") or in c("lenStruc","lenVar","ageStruc","ageVar") 
lStruc <- Slot%in%c("lenStruc","lenVar","ageStruc","ageVar")
vrbl <- ""

if (lStruc) {

vrbl <- switch(Slot,
               lenStruc="length",
               lenVar="length",
               ageStruc="age",
               ageVar="age")
#levels of 'vrbl' field must be defined as numerics
val <- as.numeric(as.character(tab[,vrbl]))
if (!is.na(step)) {
tab[,vrbl] <- factor(as.character(val),levels=seq(min(val,na.rm=TRUE),max(val,na.rm=TRUE),by=step))
} else {
tab[,vrbl] <- factor(val)
}

Xstratum <- NULL

if (is.null(dots$xlab)) dots$xlab <- vrbl 
} 


if (is.null(dots$ylab)) dots$ylab <- Slot
tstSp <- !is.na(object@species)
tstCat <- !is.na(object@catchCat)
if (is.null(dots$main)) dots$main <- paste("\"dbeOutput\" graph of '",Slot,"' slot \n",
    paste(c("for ",paste("\"",object@species,"\" species",sep="")," and ",
    paste("\"",object@catchCat,"\" fraction",sep=""))[c(tstSp|tstCat,tstSp,tstSp&tstCat,tstCat)],collapse=""),sep="")   

strip.col <- trellis.par.get("strip.background")$col

#according to 'type', the call will be quite different
plotFun <- switch(type,
                  bar="barchart",
                  line="xyplot",
                  point="xyplot")
typePar <- switch(type,
                 bar="",
                 line="type=\"l\",",
                 point="type=\"p\",") 


  #-----------------------------------------------------------------------------
  # Graphical display
  #-----------------------------------------------------------------------------

indStr <- c(timeStrata,spaceStrata,techStrata)
intLeg <- c(object@strataDesc@timeStrata,object@strataDesc@spaceStrata,object@strataDesc@techStrata)

if (all(is.na(intLeg))) intLeg <- c("time","space","technical") 

if (is.null(Xstratum)) { 

nTst <- sum(indStr)==0

eval(parse('',text=paste(plotFun,"(value ~ ",c("rep(0,nrow(tab))",vrbl)[c(!lStruc,lStruc)],paste("|",paste(c("time","space","technical")[indStr],collapse="*",sep=""),sep="")[!nTst],
  ",data=tab,horizontal=FALSE,",typePar,"drop.unused.levels=FALSE,"[lStruc & !is.na(step)],
  "main=list(dots$main,font=dots$font.main,col=dots$col.main,cex=dots$cex.main),ylim=c(0,max(tab$value,na.rm=TRUE)*1.02),",
  "col=dots$col,lwd=dots$lwd,lty=dots$lty,pch=dots$pch,cex=dots$p.cex,fill=dots$col,par.strip.text=list(font=dots$font.lab),",
  "scales=list(x=list(rot=dots$rot,cex=dots$cex",",draw=FALSE"[!lStruc],"),y=list(cex=dots$cex),font=dots$font.axis,col=dots$col.axis,cex=dots$cex.axis),",
  "xlab=list(dots$xlab,font=dots$font.lab,col=dots$col.lab,cex=dots$cex.lab),ylab=list(dots$ylab,font=dots$font.lab,col=dots$col.lab,cex=dots$cex.lab)",
  ",key=list(points=list(pch=15,cex=dots$k.cex,col=strip.col[1:sum(indStr)]),text=list(intLeg[indStr]),space=\"right\",font=dots$font.lab,columns=1,border=TRUE)"[dispKey & !nTst],
  ")",sep="")))                 

} else {

Str <- c("time","space","technical")[indStr] ; STR <- intLeg[indStr]
if (!Xstratum%in%Str) stop("Wrong 'Xstratum' parameter!!")
if (is.null(dots$xlab)) dots$xlab <- STR[match(Xstratum,Str)] 
newStr <- Str[-match(Xstratum,Str)] ; newSTR <- STR[-match(Xstratum,Str)]
dots$
indStr[match(Xstratum,c("time","space","technical"))] <- FALSE
nTst <- sum(indStr)==0

eval(parse('',text=paste(plotFun,"(value ~ ",Xstratum,paste("|",paste(newStr,collapse="*",sep=""),sep="")[!nTst],",data=tab,horizontal=FALSE,",typePar,
  "main=list(dots$main,font=dots$font.main,col=dots$col.main,cex=dots$cex.main),ylim=c(0,max(tab$value,na.rm=TRUE)*1.02),",
  "col=dots$col,lwd=dots$lwd,lty=dots$lty,pch=dots$pch,cex=dots$p.cex,fill=dots$col,par.strip.text=list(font=dots$font.lab),",
  "scales=list(x=list(rot=dots$rot,cex=dots$cex),y=list(cex=dots$cex),font=dots$font.axis,col=dots$col.axis,cex=dots$cex.axis),",
  "xlab=list(dots$xlab,font=dots$font.lab,col=dots$col.lab,cex=dots$cex.lab),ylab=list(dots$ylab,font=dots$font.lab,col=dots$col.lab,cex=dots$cex.lab)",
  ",key=list(points=list(pch=15,cex=dots$k.cex,col=strip.col[1:length(newStr)]),text=list(newSTR),space=\"right\",font=dots$font.lab,columns=1,border=TRUE)"[dispKey & !nTst],
  ")",sep=""))) 

}
})

          
###############################################################################################
###############################################################################################
###############################      EXAMPLE      #############################################
###############################################################################################
###############################################################################################
                                                                               



##usage of 'dbePlot'
#dbePlot(newObj1,"totalW",type="point",col="gold",pch=15,p.cex=1.7)
#windows()
#dbePlot(newObj1,"totalW",type="bar",Xstratum="time",cex=0.8)
#windows()
#dbePlot(newObj1,"totalW",type="bar",Xstratum="technical")
#
##creation of a 'lenStruc' slot
#norm <- dnorm(7:25,15,6)
#newObj1@lenStruc$estim <- do.call("rbind",lapply(1:19,function(x) {df <- newObj1@totalW$estim
#                                         df$length <- as.character(6+x)
#                                         df$value <- norm[x]*df$value
#                                         return(df)}))
#dbePlot(newObj1,"lenStruc",type="bar",col="gold",pch=15,p.cex=1.7)
#
#



