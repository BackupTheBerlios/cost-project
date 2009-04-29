################################################################################
# Methods to be exported
################################################################################


setGeneric("totVolume", function(dbeOutput,
                                 csObject,
                                 ceObject,
                                 clObject,
                                 ...){
	standardGeneric("totVolume")}
)



setMethod("totVolume", signature(dbeOutput="dbeOutput",csObject="csDataCons",ceObject="ceDataCons",clObject="missing"), function(dbeOutput,
                                                                                                                                 csObject,
                                                                                                                                 ceObject,
                                                                                                                                 type="trip",   #or "fo", "fd", "time"
                                                                                                                                 val="weight",  #or "number" or "nAtLength"
                                                                                                                                 sampPar=TRUE,
                                                                                                                                 ...){
if (type=="landings") stop("'landings' type requires a cs, a ce and a cl object!!")
eval(parse('',text=paste("procRaise.",type,"(csObject,ceObject,dbeOutput,val=val,sampPar=sampPar)",sep="")))

})


setMethod("totVolume", signature(dbeOutput="dbeOutput",csObject="csDataCons",ceObject="ceDataCons",clObject="clDataCons"), function(dbeOutput,
                                                                                                                                    csObject,
                                                                                                                                    ceObject,
                                                                                                                                    clObject,
                                                                                                                                    landSpp=as.character(NA),
                                                                                                                                    val="weight",  #or "number" or "nAtLength"
                                                                                                                                    sampPar=TRUE,
                                                                                                                                    ...){

procRaise.landings(csObject,ceObject,clObject,dbeOutput,landSpp=landSpp,val=val,sampPar=sampPar)

})

