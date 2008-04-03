#####################
##                 ##
##  Summary(plus)  ##
##                 ##
## MM  07/02/2008  ##
#####################



setGeneric("SummarY", function(object,
                               tab="tr",
                               sizeMax=20,
                               except=NULL,
                               ...){
	standardGeneric("SummarY")}
)


################################################################################
################################################################################


proc.SummarY <- function(dat,
                         sizeMax,
                         except=NULL,
                         indFactor,...) {   #subprocedure

#'indFactor' fields as factors
sapply(indFactor,function(x) dat[,x] <<- factor(as.character(dat[,x])))
dat <- dat[,!names(dat)%in%except]
classes <- sapply(1:ncol(dat),function(x) class(dat[,x]))
vsum1 <- summary(dat,maxsum=nrow(dat)+2) ; vsum2 <- summary(dat,maxsum=sizeMax)
index <- (1:ncol(vsum2))[classes=="factor"]

#sorting process
ord <- function(x) {val <- as.numeric(sapply(x,function(x) strsplit(x,":")[[1]][2]))
                    val[val==0] <- x[val==0] <- NA  ; indNA <- grep("NA's",x)
                    indOth <- grep("(Other)",x) ; val[indNA] <- val[indOth] <- 0
                    indic <- order(val,decreasing=TRUE) ; return(x[indic])}
                    
invisible(sapply(index,function(x) vsum1[,x] <<- ord(vsum1[,x])))
invisible(sapply(index,function(x) vsum2[,x] <<- ord(vsum2[,x])))

#take off useless lines
test1 <- apply(vsum1,1,function(x) all(is.na(x))) ; test2 <- apply(vsum2,1,function(x) all(is.na(x)))
vsum1 <- as.table(vsum1[!test1,]) ; vsum2 <- as.table(vsum2[!test2,])

#according to output table size
if (nrow(vsum1)<=sizeMax) return(vsum1) else {warning(paste("Complete summary tab size exceeds sizeMax!! (",
                                                      as.character(nrow(vsum1))," rows) ",sep=""))
                                              return(vsum2)}
}


################################################################################
################################################################################



setMethod("SummarY", signature(object="csData"), function(object,
                                                          tab="tr",
                                                          sizeMax=20,
                                                          except=NULL,...){
eval(parse('',text=paste("dat <- object@",tab,sep="")))
listFact <- list(tr=1:15,hh=c(1:14,19:29),sl=c(1:13,16),hl=1:15,ca=c(1:20,22:23,25,27:28))
proc.SummarY(dat,sizeMax,except,listFact[[tab]])                                       
})




setMethod("SummarY", signature(object="clData"), function(object,
                                                          tab="missing",
                                                          sizeMax=20,
                                                          except=NULL,...){
dat <- object@cl
proc.SummarY(dat,sizeMax,except,1:14)
})




setMethod("SummarY", signature(object="ceData"), function(object,
                                                          tab="missing",
                                                          sizeMax=20,
                                                          except=NULL,...){
dat <- object@ce
proc.SummarY(dat,sizeMax,except,1:9)
})

