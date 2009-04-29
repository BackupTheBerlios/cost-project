#-------------------------------------------------------------------------------
# costData, costDataVal, costDataCons classes and creators
# Dorleta Garcia :: Azti-Tecnalia
# 29/04/2009 15:50:50
#-------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# costData
#-------------------------------------------------------------------------------

setClass("costData",
	representation(
        desc        = "character",                     # descriptor
        ce          = "ceData",
        cl          = "clData",
		cs          = "csData"),
 	prototype(
        desc        = "costDataObject",
        ce          = ceData(),
        cl          = clData(),
		cs          = csData())
)


# missing
setGeneric("costData", function(ce,cl, cs,...){
	standardGeneric("costData")
})
setMethod("costData", signature("missing", "missing", "missing"), function(ce, cl, cs, desc="Unknown stock",
                                                                  ...){
        res      <-   new('costData')

        return(res)})

# ce, cl, cs

setMethod("costData", signature("ceData", "clData", "csData"), function(ce, cl, cs,
                                                                  desc="Unknown stock",
                                                                  ...){
        res      <- costData()
        res@desc <- desc
        res@cs   <- cs
        res@cl   <- cl
        res@ce   <- ce

        return(res)})


#--------------------------------------------------------------------------------
# costDataVal
#-------------------------------------------------------------------------------
setClass("costDataVal",
	representation(
        desc        = "character",                     # descriptor
        ce          = "ceDataVal",
        cl          = "clDataVal",
		cs          = "csDataVal"),
 	prototype(
        desc        = "costDataValObject",
        ce          = ceDataVal(),
        cl          = clDataVal(),
		cs          = csDataVal())
)

# missing
setGeneric("costDataVal", function(obj,...){
	standardGeneric("costDataVal")
})


# costData
setMethod("costDataVal", signature("costData"), function(obj, desc,
                                                                  ...){
         res     <- new('costDataVal')
         res@desc <- ifelse(missing(desc), obj@desc, desc)
         res@ce  <- ceDataVal(obj@ce)
         res@cl  <- clDataVal(obj@cl)
         res@cs  <- csDataVal(obj@cs)
         
        return(res)})


#--------------------------------------------------------------------------------
# costDataCons
#-------------------------------------------------------------------------------
setClass("costDataCons",
	representation(
        desc        = "character",                     # descriptor
        ce          = "ceDataCons",
        cl          = "clDataCons",
		cs          = "csDataCons"),
 	prototype(
        desc        = "costDataObjectCons",
        ce          = ceDataCons(),
        cl          = clDataCons(),
		cs          = csDataCons())
)

# missing
setGeneric("costDataCons", function(object, objStrat, ...){
	standardGeneric("costDataCons")
})


setMethod("costDataCons", signature("missing", 'missing'), function(desc="Unknown stock",
                                                                  ...){
        res      <- new('costDataCons')
        res@desc <- desc
        return(res)})


setMethod("costDataCons", signature("costDataVal", "strIni"), function(object,
                                                                  objStrat,
                                                                  desc="Consolidated Data",
                                                                  ...){
        res <- costDataCons()
        res@desc <- desc
        res@cs  <- csDataCons(object@cs, objStrat, desc,...)
        res@cl  <- clDataCons(object@cl, objStrat, desc,...)
        res@ce  <- ceDataCons(object@ce, objStrat, desc,...)

        return(res)})
