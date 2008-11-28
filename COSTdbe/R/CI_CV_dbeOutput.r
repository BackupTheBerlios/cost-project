

#library(COSTdbe)                          
#load("C:/CICVenvOrig.RData")              

################################################################
#                                                              #
# CI and CV calculation from 'dbeOutput' data - object updated #
#                                                              #
################################################################



#internal functions calculating CI or CV from generic dataframe, returning formatted table : 2 methods (replicates, analytical results)

  # 1. Replicates
  ######################
  
    ####################################################################
    #                                                                  #
    # ciRepFun : calculation of stratified ci from boostrap replicates #
    #                                                                  #
    ####################################################################

  #calculation from replicates                                           #object : 'dbeOutput' object 
ciRepFun <- function(object,vrbl="l",probs=c(0.025,0.975),...) {         #vrbl = "l"(length structure), "a"(age structure), "n"(total numbers) or "w"(total weights)
#'probs' parameter must be a numeric with 2 elements                     #... 'quantile' function parameters
if (!is.numeric(probs) | length(probs)!=2) stop("'probs' parameter must be a numeric with 2 elements!!")
if (probs[1]>probs[2]) stop("wrong 'probs' parameter!!")

#relevant data according to 'vrbl' parameter
  #'rep' data
dfRep <- switch(vrbl,
                l=object@lenStruc$rep,
                a=object@ageStruc$rep,
                n=object@totalN$rep,
                w=object@totalW$rep)

  #is 'vrbl' correctly defined ?
if (is.null(dfRep)) stop("wrong 'vrbl' parameter!!")
if (all(is.na(dfRep))) stop("missing data in input object!!")
  #is there any data?
if (nrow(dfRep)==0) stop("no available data!!")

#definition of aggregation fields
exprBy <- paste("length=as.numeric(as.character(dfRep$length)),"["l"%in%vrbl],
                "age=as.numeric(as.character(dfRep$age)),"["a"%in%vrbl],
                "technical=dfRep$technical,"[!all(is.na(dfRep$technical))],
                "space=dfRep$space,"[!all(is.na(dfRep$space))],
                "time=dfRep$time,"[!all(is.na(dfRep$time))],
                sep="",collapse="")
#last character is removed
exprBy <- substr(exprBy,1,nchar(exprBy)-1)
#stratified CI
CIdf1 <- CIdf2 <- NULL
eval(parse('',text=paste("CIdf1 <- aggregate(dfRep$value,list(",exprBy,"),function(x) quantile(x,probs=probs,...)[1])",sep=""))) ; names(CIdf1)[ncol(CIdf1)] <- "inf" 
eval(parse('',text=paste("CIdf2 <- aggregate(dfRep$value,list(",exprBy,"),function(x) quantile(x,probs=probs,...)[2])",sep=""))) ; names(CIdf2)[ncol(CIdf2)] <- "sup" 
eval(parse('',text=paste("CIdf3 <- aggregate(dfRep$value,list(",exprBy,"),mean,na.rm=TRUE)",sep=""))) ; names(CIdf3)[ncol(CIdf3)] <- "value" 
#result
df <- merge(CIdf1,CIdf2,sort=FALSE)
df <- merge(df,CIdf3,sort=FALSE)

#formatting process
DF <- dfRep[rep(1,nrow(df)),-match(c("iter"),names(dfRep))]
DF$sup <- DF$inf <- DF$value <- 0
invisible(sapply(names(df),function(x) DF[,x] <<- df[,x]))
rownames(DF) <- 1:nrow(DF)

return(DF)
}




    ####################################################################
    #                                                                  #
    # cvRepFun : calculation of stratified cv from boostrap replicates #
    #                                                                  #
    ####################################################################

  #calculation from replicates                      #object : 'dbeOutput' object                      
cvRepFun <- function(object,vrbl="l",...) {         #vrbl = "l"(length structure), "a"(age structure), "n"(total numbers) or "w"(total weights)
                                                    #... 'quantile' function parameters
#relevant data according to 'vrbl' parameter
  #'rep' data
dfRep <- switch(vrbl,
                l=object@lenStruc$rep,
                a=object@ageStruc$rep,
                n=object@totalN$rep,
                w=object@totalW$rep)

  #is 'vrbl' correctly defined ?
if (is.null(dfRep)) stop("wrong 'vrbl' parameter!!")
if (all(is.na(dfRep))) stop("missing data in input object!!")
  #is there any data?
if (nrow(dfRep)==0) stop("no available data!!")

#definition of aggregation fields
exprBy <- paste("length=as.numeric(as.character(dfRep$length)),"["l"%in%vrbl],
                "age=as.numeric(as.character(dfRep$age)),"["a"%in%vrbl],
                "technical=dfRep$technical,"[!all(is.na(dfRep$technical))],
                "space=dfRep$space,"[!all(is.na(dfRep$space))],
                "time=dfRep$time,"[!all(is.na(dfRep$time))],
                sep="",collapse="")
#last character is removed
exprBy <- substr(exprBy,1,nchar(exprBy)-1)
#stratified CV
CVdf <- NULL
eval(parse('',text=paste("CVdf <- aggregate(dfRep$value,list(",exprBy,"),function(x) sd(x)/mean(x))",sep=""))) 
names(CVdf)[ncol(CVdf)] <- "value" 

#formatting process
DF <- dfRep[rep(1,nrow(CVdf)),-match(c("iter"),names(dfRep))]
invisible(sapply(names(CVdf),function(x) DF[,x] <<- CVdf[,x]))
rownames(DF) <- 1:nrow(DF)

return(DF)
}





  # 2. Estimates
  ######################



    ##########################################################
    #                                                        #
    # ciEstFun : calculation of stratified ci from estimates #
    #                                                        #
    ##########################################################

  #calculation from estimates                                            #object : 'dbeOutput' object 
ciEstFun <- function(object,vrbl="l",probs=c(0.025,0.975),...) {         #vrbl = "l"(length structure), "a"(age structure), "n"(total numbers) or "w"(total weights)
#'probs' parameter must be a numeric with 2 elements
if (!is.numeric(probs) | length(probs)!=2) stop("'probs' parameter must be a numeric with 2 elements!!")
if (probs[1]>probs[2]) stop("wrong 'probs' parameter!!")

#relevant data according to 'vrbl' parameter
  #"estim" data
dfEstim <- switch(vrbl,
                l=object@lenStruc$estim,
                a=object@ageStruc$estim,
                n=object@totalN$estim,
                w=object@totalW$estim)
               
  #'var' data
dfVar <- switch(vrbl,
                l=object@lenVar,
                a=object@ageVar,
                n=object@totalNvar,
                w=object@totalWvar)

  #is 'vrbl' correctly defined ?
if (is.null(dfEstim) | is.null(dfVar)) stop("wrong 'vrbl' parameter!!")
if (all(is.na(dfEstim)) | all(is.na(dfVar))) stop("missing data in input object!!")
  #is there any data?
if (nrow(dfEstim)==0 | nrow(dfVar)==0) stop("missing data in input object!!")
  #are the tables consistent?
if (!identical(dim(dfEstim),dim(dfVar))) stop("tables are not matching!!") 

#both tables are merged
names(dfVar)[ncol(dfVar)] <- "var"
df <- merge(dfEstim,dfVar)
  #is number of rows the same in df than in previous tables?
if (!identical(nrow(dfEstim),nrow(df))) stop("tables are not matching!!") 

#'inf' and 'sup' fields are calulated, and pasted to df
qInd <- qnorm(probs)
dfInfSup <- data.frame(df$value+t(qInd%*%t(sqrt(df$var)))) ; names(dfInfSup) <- c("inf","sup")
df <- cbind(df,dfInfSup)

#problems of negative bounds
if (any(c(df$inf,df$sup)<0)) warning("negative CI bound(s)!!")                  

#formatting process
DF <- df[,c(names(dfEstim),"inf","sup")]
rownames(DF) <- 1:nrow(DF)

return(DF)
}




    ##########################################################
    #                                                        #
    # cvEstFun : calculation of stratified cv from estimates #
    #                                                        #
    ##########################################################

  #calculation from estimates                       #object : 'dbeOutput' object                    
cvEstFun <- function(object,vrbl="l",...) {         #vrbl = "l"(length structure), "a"(age structure), "n"(total numbers) or "w"(total weights)

#relevant data according to 'vrbl' parameter
  #"estim" data
dfEstim <- switch(vrbl,
                l=object@lenStruc$estim,
                a=object@ageStruc$estim,
                n=object@totalN$estim,
                w=object@totalW$estim)
                
  #'var' data
dfVar <- switch(vrbl,
                l=object@lenVar,
                a=object@ageVar,
                n=object@totalNvar,
                w=object@totalWvar)
                            
  #is 'vrbl' correctly defined ?
if (is.null(dfEstim) | is.null(dfVar)) stop("wrong 'vrbl' parameter!!")
if (all(is.na(dfEstim)) | all(is.na(dfVar))) stop("missing data in input object!!")
  #is there any data?
if (nrow(dfEstim)==0 | nrow(dfVar)==0) stop("missing data in input object!!")
  #are the tables consistent?
if (!identical(dim(dfEstim),dim(dfVar))) stop("tables are not matching!!") 

#both tables are merged
nam <- names(dfEstim)
names(dfEstim)[ncol(dfEstim)] <- "estim"
names(dfVar)[ncol(dfVar)] <- "var"
df <- merge(dfEstim,dfVar)
  #is number of rows the same in df than in previous tables?
if (!identical(nrow(dfEstim),nrow(df))) stop("tables are not matching!!") 

#CV is calulated, and inserted in df
df$value <- sqrt(df$var)/df$estim

#formatting process
DF <- df[,nam]
rownames(DF) <- 1:nrow(DF)

return(DF)
}


#####################################################################################
#####################################################################################
#####################                                         #######################
#####################                 Methods                 #######################
#####################              dbeCI - dbeCV              #######################
#####################                                         #######################
#####################################################################################
#####################################################################################




        
setGeneric("dbeCalc", function(object,              # 'dbeOutput' object
                               type="CI",           # "CI" for confidence interval, or "CV" for coefficient of variation
                               vrbl="l",            # specifies data on which calculation is applied : "l" for length structure, "a" for age structure, "n" for total number estimates, "w" for total weight estimates
                               probs=c(0.025,0.975),# used only if type="CI", defines bounds
                               replicates=FALSE,    # if TRUE, calculation is made from $rep elements ; if FALSE, $estim and @...Var are used 
                               update=FALSE,        # if TRUE, updated 'dbeOutput' object is returned ; if FALSE, only resulting dataframe is returned
                               ...){
standardGeneric("dbeCalc")
})



 
setMethod("dbeCalc",signature(object="dbeOutput"),function(object,              #'dbeOutput' object 
                                                           type="CI",           # "CI" for confidence interval, or "CV" for coefficient of variation
                                                           vrbl="l",            # specifies data on which calculation is applied : "l" for length structure, "a" for age structure, "n" for total number estimates, "w" for total weight estimates
                                                           probs=c(0.025,0.975),# used only if type="CI", defines bounds
                                                           replicates=FALSE,    # if TRUE, calculation is made from $rep elements ; if FALSE, $estim and @...Var are used 
                                                           update=FALSE,        # if TRUE, updated 'dbeOutput' object is returned ; if FALSE, only resulting dataframe is returned
                                                           ...){                                

output <- NULL
#according to input parameters, one of the previous functions is to be used
lType <- tolower(type)
if (!lType%in%c("ci","cv")) stop("wrong 'type' parameter!!") 
if (!is.logical(replicates)) stop("wrong 'replicates' parameter!!") 
if (!is.logical(update)) stop("wrong 'update' parameter!!")

if (replicates) lData <- "Rep" else lData <- "Est" 

#so, expected table is...
eval(parse('',text=paste("output <- ",lType,lData,"Fun(object=object,vrbl=vrbl,probs=probs,...)",sep=""))) 

#output depends on 'update' parameter
if (update) {
  #slot to update must be identified 
  updSlot <- switch(vrbl,
                    l="lenNum",
                    a="ageNum",
                    n="totalNnum",
                    w="totalWnum")
  #and then, object is updated and returned
  eval(parse('',text=paste("object@",updSlot,"$",lType," <- output",sep=""))) 
  return(object)
} else {
  return(output)
}

})


#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

############ 
# Examples #
############ 
 
#dbeCalc(object) 
#dbeCalc(object,type="CV") 
#dbeCalc(object,type="CI",vrbl="a") 
#dbeCalc(object,type="CI",vrbl="n",probs=c(0.1,0.9)) #80% 
#object <- dbeCalc(object,type="ci",replicates=TRUE,update=TRUE) 
#object <- dbeCalc(object,type="cv",vrbl="a",replicates=TRUE,update=TRUE) 
#object <- dbeCalc(object,type="ci",vrbl="n",update=TRUE) 
#object <- dbeCalc(object,type="cv",vrbl="n",update=TRUE) 
 
 
 
 
 
                                                                   