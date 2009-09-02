#-------------------------------------------------------------------------------
# Function to combine output objects  
# Dorleta Garcia - Azti-Tecnalia
# 28/08/2009 12:23:47
#-------------------------------------------------------------------------------

# objs: a list with name of the objects to merge.
# name: file root name
# seq:file index
# nsample: number of simulated samples in each file.
# dir: Full path of the directory where the files are located, if NULL the 
#       files are loaded from the current R directory. 

combOutSim <- function(objs, name, seq, nsample, dir = NULL, rep = FALSE){
    
        if(!is.null(dir))  setwd(dir)
        
        slots <- list(lenStruc = c('estim', 'rep'), lenVar = NA,
                      ageStruc = c('estim', 'rep'), ageVar = NA,
                      totalN   = c('estim', 'rep'), totalN = NA,
                      totalW   = c('estim', 'rep'), totalW = NA)
        
        res        <- vector('list', length(objs))
        names(res) <- names(objs)
        
        load(paste(seq[1], name, sep = "_"))
        
        for(j in 1:length(objs))  res[[j]] <- get(objs[j])
        
        k <- 1
        for(i in seq[-1]){
            cat('-------------------------------------------------------------\n')
            cat('File: ', i, '\n')
            
            load(paste(i, name, sep = "_"))
            
            for(j in 1:length(objs)){
            
                for(s in 1:length(slots)){

                    x <- slot(res[[j]], names(slots)[s])
                    
                    if(class(x) != 'list'){
                        y        <- slot(get(objs[j]), names(slots)[s])
                        y$sample <- as.character(k*2 + as.numeric(y$sample))
                        slot(res[[j]], names(slots)[s]) <- rbind(slot(res[[j]], names(slots)[s]), y)
                    }
                    else{
                        y        <- slot(get(objs[j]), names(slots)[s])$estim
                        y$sample <- as.character(k*2 + as.numeric(y$sample))
                        slot(res[[j]], names(slots)[s])$estim <- rbind(slot(res[[j]], names(slots)[s])$estim, y)
                        
                        if(rep == TRUE){
                            y        <- slot(get(objs[j]), names(slots)[s])$rep
                            y$sample <- as.character(k*2 + as.numeric(y$sample))
                            slot(res[[j]], names(slots)[s])$rep <- rbind(slot(res[[j]], names(slots)[s])$rep, y)
                       }
                }
             }
              
        }
         k <- k+1
    }
    return(res)
}           

