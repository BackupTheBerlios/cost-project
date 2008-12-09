#####################
##                 ##
##  Summary(plus)  ##
##                 ##
## MM  07/02/2008  ##
#####################



setGeneric("cstSummary", function(object,
                                  tab="tr",
                                  sizeMax=20,
                                  except=NULL,
                                  ...){
	standardGeneric("cstSummary")}
)


################################################################################
################################################################################


proc.cstSummary <- function(dat,
                            sizeMax,
                            except=NULL,
                            indFactor,...) {   #subprocedure

#'indFactor' fields as factors
sapply(indFactor,function(x) dat[,x] <<- factor(as.character(dat[,x])))
dat <- dat[,!names(dat)%in%except]
classes <- sapply(1:ncol(dat),function(x) class(dat[,x]))
vsum1 <- summary(dat,maxsum=nrow(dat)+2)
vsum2 <- summary(dat,maxsum=sizeMax)
index <- (1:ncol(vsum2))[classes=="factor"]

#sorting process
ord <- function(x) {val <- as.numeric(sapply(x,function(x) strsplit(x,":")[[1]][2]))
                    val[val==0] <- x[val==0] <- NA
                    indNA <- grep("NA's",x)
                    indOth <- grep("(Other)",x)
                    val[indNA] <- val[indOth] <- 0
                    indic <- order(val,decreasing=TRUE)
                    return(x[indic])}
                    
invisible(sapply(index,function(x) vsum1[,x] <<- ord(vsum1[,x])))
invisible(sapply(index,function(x) vsum2[,x] <<- ord(vsum2[,x])))

#take off useless lines
test1 <- apply(vsum1,1,function(x) all(is.na(x)))
test2 <- apply(vsum2,1,function(x) all(is.na(x)))
vsum1 <- as.table(vsum1[!test1,])
vsum2 <- as.table(vsum2[!test2,])

#according to output table size
if (nrow(vsum1)<=sizeMax) 
  return(vsum1) 
else {
  warning(paste("Complete summary tab size exceeds sizeMax!! (",as.character(nrow(vsum1))," rows) ",sep=""))
  return(vsum2)}
}


################################################################################
################################################################################



setMethod("cstSummary", signature(object="csData"), function(object,
                                                             tab="tr",
                                                             sizeMax=20,
                                                             except=NULL,
                                                             ...){
eval(parse('',text=paste("dat <- object@",tab,sep="")))
listFact <- list(tr=1:16,hh=c(1:14,19:29),sl=c(1:14,17),hl=1:15,ca=c(1:21,23:25,27,29:31))
proc.cstSummary(dat,sizeMax,except,listFact[[tab]])                                       
})



setMethod("cstSummary", signature(object="csDataVal"), function(object,
                                                                tab="tr",
                                                                sizeMax=20,
                                                                except=NULL,
                                                                ...){
eval(parse('',text=paste("dat <- object@",tab,sep="")))
listFact <- list(tr=1:16,hh=c(1:14,19:29),sl=c(1:14,17),hl=1:15,ca=c(1:21,23:25,27,29:31))
proc.cstSummary(dat,sizeMax,except,listFact[[tab]])                                       
})



setMethod("cstSummary", signature(object="csDataCons"), function(object,
                                                                 tab="tr",
                                                                 sizeMax=20,
                                                                 except=NULL,
                                                                 ...){
eval(parse('',text=paste("dat <- object@",tab,sep="")))
listFact <- list(tr=2:14,hh=3:16,sl=c(4:15,18),hl=4:16,ca=c(3:17,19:21,23,25:27))
proc.cstSummary(dat,sizeMax,except,listFact[[tab]])                                       
})






setMethod("cstSummary", signature(object="clData"), function(object,
                                                             tab="missing",
                                                             sizeMax=20,
                                                             except=NULL,
                                                             ...){
dat <- object@cl
proc.cstSummary(dat,sizeMax,except,1:17)
})



setMethod("cstSummary", signature(object="clDataVal"), function(object,
                                                                tab="missing",
                                                                sizeMax=20,
                                                                except=NULL,
                                                                ...){
dat <- object@cl
proc.cstSummary(dat,sizeMax,except,1:17)
})



setMethod("cstSummary", signature(object="clDataCons"), function(object,
                                                                 tab="missing",
                                                                 sizeMax=20,
                                                                 except=NULL,...){
dat <- object@cl
proc.cstSummary(dat,sizeMax,except,1:9)
})





setMethod("cstSummary", signature(object="ceData"), function(object,
                                                             tab="missing",
                                                             sizeMax=20,
                                                             except=NULL,
                                                             ...){
dat <- object@ce
proc.cstSummary(dat,sizeMax,except,1:12)
})



setMethod("cstSummary", signature(object="ceDataVal"), function(object,
                                                                tab="missing",
                                                                sizeMax=20,
                                                                except=NULL,
                                                                ...){
dat <- object@ce
proc.cstSummary(dat,sizeMax,except,1:12)
})



setMethod("cstSummary", signature(object="ceDataCons"), function(object,
                                                                 tab="missing",
                                                                 sizeMax=20,
                                                                 except=NULL,
                                                                 ...){
dat <- object@ce
proc.cstSummary(dat,sizeMax,except,1:4)
})




################################################################################################################################
# tapply-like methods to display information from csDataVal, clDataVal, ceDataVal, csDataCons, clDataCons & ceDataCons objects #
################################################################################################################################


setGeneric("disInfo", function(object,...){
	standardGeneric("disInfo")
	}
)


################################
################################
##  Validated data structure  ##
################################
################################



setMethod("disInfo", signature("csDataVal"), function(object,path,field,by,fun,...,biopar=FALSE,title="",append=TRUE) {

if (missing(path)) stop("Missing 'path' parameter!!")
if (missing(field)) stop("Missing 'field' parameter!!")
if (length(field)==0) stop("Missing 'field' parameter!!")
if (missing(by)) stop("Missing 'by' parameter!!")
if (missing(fun)) stop("Missing 'fun' parameter!!")

allFields <- c(field,by)
if (biopar) {
  TAB <- ca(object)
  if (any(!allFields%in%names(TAB))) stop("some parameters don't match with CA table!!") 
} else {
#"month", "quarter" & "semester" are put in HH
  HH <- hh(object)
  HH$month <- as.numeric(sapply(HH$date,function(x) strsplit(as.character(x),"-")[[1]][2]))
  HH$quarter <- ceiling(HH$month/3) ; HH$semester <- ceiling(HH$month/6)
#step by step, we merge all the tables from tr to hl until all specified fields are in
  go <- TRUE
#1st step
  tab <- tr(object)
  if (all(allFields%in%names(tab))) go <- FALSE
#2nd step
  if (go) {
    tab <- tab[tab$sampType!="V",]
    tab <- merge(HH,tab,by=c("sampType","landCtry","vslFlgCtry","year","proj","trpCode"),sort=FALSE,all=TRUE)
  }
  if (all(allFields%in%names(tab))) go <- FALSE
#3rd step
  if (go) tab <- merge(sl(object),tab,by=c("sampType","landCtry","vslFlgCtry","year","proj","trpCode","staNum"),sort=FALSE,all.x=TRUE)
  if (all(allFields%in%names(tab))) go <- FALSE
#4th step
  if (go) tab <- merge(slSex(sl(object),hl(object)),tab,by=c("sampType","landCtry","vslFlgCtry","year","proj","trpCode",
                                           "staNum","spp","catchCat","landCat","commCatScl","commCat","subSampCat","sex"),sort=FALSE,all.x=TRUE)
  TAB <- tab
  if (any(!allFields%in%names(TAB))) stop("some parameters don't match with tables!!") 
  }

if (length(field)==1) {
  eval(parse('',text=paste("calc <- tapply(TAB$",field,",list(",paste("TAB$",by,collapse=",",sep=""),"),fun,...)",sep="")))
} else {
  Field <- apply(TAB[,field],1,paste,collapse=":-:")
  eval(parse('',text=paste("calc <- tapply(Field,list(",paste("TAB$",by,collapse=",",sep=""),"),fun,...)",sep="")))
}

#Report
sink(file=path,append=append)
    cat("\n") 		
    cat(title,"\n")
    eval(parse('',text="cat(paste(rep(\"-\",nchar(title)),collapse=\"\",sep=\"\"),\"\n\")"[nchar(title)>0]))
		cat("\n")
    cat("Field -> ",paste(field,collapse=", ",sep=""),"\n")
		cat("By    -> ",paste(by,collapse=", ",sep=""),"\n")
		cat("FUN   -> ",deparse(substitute(fun)),"\n","\n")
    print(calc)
		sink()
invisible(calc)
})






setMethod("disInfo", signature("clDataVal"), function(object,path,field,by,fun,...,title="",append=TRUE) {

if (missing(path)) stop("Missing 'path' parameter!!")
if (missing(field)) stop("Missing 'field' parameter!!")
if (length(field)!=1) stop("One field must be specified in 'field' parameter!!")
if (missing(by)) stop("Missing 'by' parameter!!")
if (missing(fun)) stop("Missing 'fun' parameter!!")

tab <- cl(object)
allFields <- c(field,by)
if (!all(allFields%in%names(tab))) stop("Fields don't match with 'cl' format!!")
eval(parse('',text=paste("calc <- tapply(tab$",field,",list(",paste("tab$",by,collapse=",",sep=""),"),fun,...)",sep="")))
#Report
sink(file=path,append=append)
    cat("\n") 		
    cat(title,"\n")
    eval(parse('',text="cat(paste(rep(\"-\",nchar(title)),collapse=\"\",sep=\"\"),\"\n\")"[nchar(title)>0]))
		cat("\n")
    cat("Field -> ",field,"\n")
		cat("By    -> ",paste(by,collapse=", ",sep=""),"\n")
		cat("FUN   -> ",deparse(substitute(fun)),"\n","\n")
    print(calc)
		sink()
invisible(calc)
})







setMethod("disInfo", signature("ceDataVal"), function(object,path,field,by,fun,...,title="",append=TRUE) {

if (missing(path)) stop("Missing 'path' parameter!!")
if (missing(field)) stop("Missing 'field' parameter!!")
if (length(field)!=1) stop("One field must be specified in 'field' parameter!!")
if (missing(by)) stop("Missing 'by' parameter!!")
if (missing(fun)) stop("Missing 'fun' parameter!!")

tab <- ce(object)
allFields <- c(field,by)
if (!all(allFields%in%names(tab))) stop("Fields don't match with 'ce' format!!")
eval(parse('',text=paste("calc <- tapply(tab$",field,",list(",paste("tab$",by,collapse=",",sep=""),"),fun,...)",sep="")))
#Report
sink(file=path,append=append)
    cat("\n") 		
    cat(title,"\n")
    eval(parse('',text="cat(paste(rep(\"-\",nchar(title)),collapse=\"\",sep=\"\"),\"\n\")"[nchar(title)>0]))
		cat("\n")
    cat("Field -> ",field,"\n")
		cat("By    -> ",paste(by,collapse=", ",sep=""),"\n")
		cat("FUN   -> ",deparse(substitute(fun)),"\n","\n")
    print(calc)
		sink()
invisible(calc)
})









###################################
###################################
##  Consolidated data structure  ##
###################################
###################################



setMethod("disInfo", signature("csDataCons"), function(object,path,field,by,fun,...,biopar=FALSE,title="",append=TRUE) {   #by only contains "time","space" or/and "technical

if (missing(path)) stop("Missing 'path' parameter!!")
if (missing(field)) stop("Missing 'field' parameter!!")
if (length(field)==0) stop("Missing 'field' parameter!!")
if (missing(by)) stop("Missing 'by' parameter!!")
by <- by[by%in%c("time","space","technical")]
if (length(by)==0) stop("'by' must contain \"time\",\"space\" or/and \"technical\"!!")
if (missing(fun)) stop("Missing 'fun' parameter!!")

if (biopar) {
  TAB <- ca(object)
  if (any(!field%in%names(TAB))) stop("'field' parameter doesn't match with CA table!!") 
} else {
  TAB <- tr(object)
  if (any(!field%in%names(TAB))) {
    TAB <- hh(object)
    if (any(!field%in%names(TAB))) {
      TAB <- sl(object)
      if (any(!field%in%names(TAB))) {
        TAB <- hl(object)
        if (any(!field%in%names(TAB))) stop("'field' parameter doesn't match with tables!!")
      }
    }
  }
}

if (length(field)==1) {
  eval(parse('',text=paste("calc <- tapply(TAB$",field,",list(",paste("TAB$",by,collapse=",",sep=""),"),fun,...)",sep="")))
} else {
  Field <- apply(TAB[,field],1,paste,collapse=":-:")
  eval(parse('',text=paste("calc <- tapply(Field,list(",paste("TAB$",by,collapse=",",sep=""),"),fun,...)",sep="")))
}

#Report
sink(file=path,append=append)
    cat("\n") 		
    cat(title,"\n")
    eval(parse('',text="cat(paste(rep(\"-\",nchar(title)),collapse=\"\",sep=\"\"),\"\n\")"[nchar(title)>0]))
		cat("\n")
    cat("Field -> ",paste(field,collapse=", ",sep=""),"\n")
		cat("By    -> ",paste(by,collapse=", ",sep=""),"\n")
		cat("FUN   -> ",deparse(substitute(fun)),"\n","\n")
    print(calc)
		sink()
invisible(calc)
})







setMethod("disInfo", signature("clDataCons"), function(object,path,field,by,fun,...,title="",append=TRUE) {   #by only contains "time","space" or/and "technical

if (missing(path)) stop("Missing 'path' parameter!!")
if (missing(field)) stop("Missing 'field' parameter!!")
if (length(field)==0) stop("Missing 'field' parameter!!")
if (length(field)>1) stop("Only one field name in 'field' parameter!!")
if (missing(by)) stop("Missing 'by' parameter!!")
by <- by[by%in%c("time","space","technical")]
if (length(by)==0) stop("'by' must contain \"time\",\"space\" or/and \"technical\"!!")
if (missing(fun)) stop("Missing 'fun' parameter!!")


TAB <- cl(object)
if (!field%in%names(TAB)) stop("'field' parameter doesn't match with CL table!!") 

eval(parse('',text=paste("calc <- tapply(TAB$",field,",list(",paste("TAB$",by,collapse=",",sep=""),"),fun,...)",sep="")))
#Report
sink(file=path,append=append)
    cat("\n") 		
    cat(title,"\n")
    eval(parse('',text="cat(paste(rep(\"-\",nchar(title)),collapse=\"\",sep=\"\"),\"\n\")"[nchar(title)>0]))
		cat("\n")
    cat("Field -> ",field,"\n")
		cat("By    -> ",paste(by,collapse=", ",sep=""),"\n")
		cat("FUN   -> ",deparse(substitute(fun)),"\n","\n")
    print(calc)
		sink()
invisible(calc)
})





setMethod("disInfo", signature("ceDataCons"), function(object,path,field,by,fun,...,title="",append=TRUE) {   #by only contains "time","space" or/and "technical

if (missing(path)) stop("Missing 'path' parameter!!")
if (missing(field)) stop("Missing 'field' parameter!!")
if (length(field)==0) stop("Missing 'field' parameter!!")
if (length(field)>1) stop("Only one field name in 'field' parameter!!")
if (missing(by)) stop("Missing 'by' parameter!!")
by <- by[by%in%c("time","space","technical")]
if (length(by)==0) stop("'by' must contain \"time\",\"space\" or/and \"technical\"!!")
if (missing(fun)) stop("Missing 'fun' parameter!!")


TAB <- ce(object)
if (!field%in%names(TAB)) stop("'field' parameter doesn't match with CE table!!") 

eval(parse('',text=paste("calc <- tapply(TAB$",field,",list(",paste("TAB$",by,collapse=",",sep=""),"),fun,...)",sep="")))
#Report
sink(file=path,append=append)
    cat("\n") 		
    cat(title,"\n")
    eval(parse('',text="cat(paste(rep(\"-\",nchar(title)),collapse=\"\",sep=\"\"),\"\n\")"[nchar(title)>0]))
		cat("\n")
    cat("Field -> ",field,"\n")
		cat("By    -> ",paste(by,collapse=", ",sep=""),"\n")
		cat("FUN   -> ",deparse(substitute(fun)),"\n","\n")
    print(calc)
		sink()
invisible(calc)
})



# CS/CL/CE consistency (presence of values of a given field in each specified table)

tabConsist <- function(lTab,field) {                 #lTab is a list containing various COST objects

newLtab <- NULL

#test on objects classes
ind <- FALSE
clas <- unlist(lapply(lTab,class))
if (all(clas%in%c("csData","ceData","clData"))) ind <- TRUE
if (all(clas%in%c("csDataVal","ceDataVal","clDataVal"))) ind <- TRUE
if (all(clas%in%c("csDataCons","ceDataCons","clDataCons"))) ind <- TRUE
if (!ind) stop("objects classes don't match!!")

#names for each table
nam <- rep(c("CS","CL","CE"),each=3) ; reit <- rep(c(5,1,1),each=3)
names(reit) <- names(nam) <- c("csData","csDataVal","csDataCons","clData","clDataVal","clDataCons","ceData","ceDataVal","ceDataCons")
NAM <- nam[clas] ; IND <- rep(1,length(NAM))
if (length(NAM)>1) invisible(sapply(2:length(NAM),function(x) IND[x] <<- sum(NAM[1:x]%in%NAM[x])))
NAM <- paste(NAM,IND,sep="")
namTab <- namT <- unlist(do.call("list",lapply(lTab,function(x) slotNames(x)[-1])))
at <- rep("_",length(namTab))
at[namTab%in%c("cl","ce")] <- ""
namTab[namTab%in%c("cl","ce")] <- ""
NAMES <- paste(rep(NAM,reit[clas]),at,namTab,sep="")

#tables in lTab are split and put in a list
eval(parse('',text=paste("newLtab <- list(",paste(NAMES,"=lTab[[",rep(1:length(NAM),reit[clas]),"]]@",namT,collapse=",",sep=""),")",sep="")))

#test the presence of specified field in the tables, and remove the table if no
n <- length(newLtab)
invisible(sapply(n:1,function(x) if (!field%in%names(newLtab[[x]])) newLtab[[x]] <<- NULL))
if (length(newLtab)==0) stop("specified field cannot be found in input objects!!")

#all values of specified field (colNames of the resulting table)
lValues <- do.call("list",lapply(newLtab,function(x) unique(x[,field])))
lValues <- sort(na.omit(unique(unlist(lValues))))

#presence of each value in each table for the specified field (Y/N)
tab <- do.call("rbind",lapply(1:length(newLtab),function(x) lValues%in%newLtab[[x]][,field]))
dimnames(tab) <- list(names(newLtab),lValues)
tab[tab] <- "x" ; tab[tab=="FALSE"] <- ""

return(as.table(tab))
}





#mix between 'aggregate' & 'tapply'
dfApply <- function(tab,valueField,rowFields,colField,fun,...){
if (length(colField)!=1) stop("wrong 'rowField' parameter!!")
if (length(valueField)!=1) stop("wrong 'valueField' parameter!!")

#rowFields are concatenated and then tapply function is used
mat <- tapply(tab[,valueField],list(apply(tab[,rowFields,drop=FALSE],1,paste,collapse=":-:"),tab[,colField]),fun,...)
mat1 <- do.call("rbind",lapply(rownames(mat),function(x) strsplit(x,":-:")[[1]]))
df1 <- as.data.frame(mat1) ; names(df1) <- rowFields
df2 <- as.data.frame(mat) ; rownames(df2) <- NULL
dfInt <- data.frame(rep("|",nrow(df1))) ; names(dfInt)[1] <- "|"
df <- cbind.data.frame(df1,dfInt,df2)
return(df)
}


