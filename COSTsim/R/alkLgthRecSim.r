#-------------------------------------------------------------------------------
# alkLgthRecSim method.
# Dorleta Garcia :: Azti-Tecnalia
# 29/04/2009 15:50:50
#-------------------------------------------------------------------------------
  
setGeneric("alkLgthRecSim", function(object,
                                  type="stepIncr",        
                                  value,                  
                                  preview=FALSE,          
                                  postview=TRUE,
                                  update=FALSE,
                                  ...) {       
	standardGeneric("alkLgthRecSim")}
)



setMethod("alkLgthRecSim",signature(object="simDataCons"), function(object,
                                                                type="stepIncr",        
                                                                value,                  
                                                                preview=FALSE,          
                                                                postview=FALSE,
                                                                update=FALSE,
                                                                ...) {
    nsamples <- length(object@samples)
    
    
    if(update == TRUE){ 
        res <- object
        for(i in 1:nsamples){
            res@samples[[i]] <- alkLgthRec(object@samples[[i]], type = type, value = value, 
                                    preview = preview, postview = postview, update = update,...)    
        }
    }
    else{
        res <- vector('list', nsamples)
         for(i in 1:nsamples){
            res[[i]]<- alkLgthRec(object@samples[[i]], type = type, value = value, 
                                    preview = preview, postview = postview, update = update,...)    
        }
    }
    
        return(res)                                                            
                                                                 
})

                 

