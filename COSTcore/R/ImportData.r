

#COSTcore package
#Read .csv files in a folder, and build consequently COST objects

importValidity <- function(df){                                       #tests on data.frames
tabnames <- c("TR","HH","SL","HL","CA","CL","CE")
RecordType <- unique(as.character(df[,1]))
if (length(RecordType)>1) stop("RecordType is not unique!!")
if (!(RecordType%in%tabnames)) stop("RecordType is not valid!!")
if (dim(df)[2]<2) stop("Wrong df size!!")
df <- df[,2:dim(df)[2]]
return(df)
}


ImportDataCost <- function(file,list.data){
df <- read.table(file,header=TRUE,sep=";")  
val <- as.character(df[1,1])
df <- importValidity(df)
eval(parse('',text=paste("list.data$",tolower(val)," <- df",sep="")))
return(list.data)
}

####################
## Main procedure ##
####################

CreateClassCOST <- function (directory,desc="My Stock") {    #directory = .csv files ; desc = slot 'desc' 
require(COSTcore)
obj <- list(cs=new("csData"),cl=new("clData"),ce=new("ceData"))
list.data <- list(tr=obj$cs@tr,hh=obj$cs@hh,sl=obj$cs@sl,hl=obj$cs@hl,ca=obj$cs@ca,cl=obj$cl@cl,ce=obj$ce@ce)
list.nms <- list(tr=names(obj$cs@tr),hh=names(obj$cs@hh),sl=names(obj$cs@sl),hl=names(obj$cs@hl),ca=names(obj$cs@ca),cl=names(obj$cl@cl),ce=names(obj$ce@ce))
Files <- paste(directory,list.files(directory),sep="/")
for (i in 1:length(Files)) list.data <- ImportDataCost(Files[i],list.data)
for (i in 1:7) names(list.data[[i]]) <- list.nms[[i]]
obj <- list(cs=new("csData",tr=list.data$tr,hh=list.data$hh,sl=list.data$sl,hl=list.data$hl,ca=list.data$ca),cl=new("clData",cl=list.data$cl),ce=new("ceData",ce=list.data$ce))
neodesc <- rep(desc,3)
obj$cs@desc <- neodesc[1] ; obj$cl@desc <- neodesc[2] ; obj$ce@desc <- neodesc[3]
invisible(obj)
}

#Example
#Put the .csv files in one folder 
directory <- "Q:/FF5_sole"   #path of the folder
obj <- CreateClassCOST(directory)
#obj$cs, obj$cl, obj$ce