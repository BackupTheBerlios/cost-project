fillGaps <- function(dbeSimObj, ageMin, ageMax, lenMin, lenMax){

    if(dbeSimObj@catchCat == 'LAN')
          res <- fillGaps.lan(dbeSimObj, ageMin, ageMax, lenMin, lenMax)
    else
        res <- fillGaps.dis(dbeSimObj, ageMin, ageMax, lenMin, lenMax)

    return(res)
}
        

fillGaps.lan <- function(dbeSimObj, ageMin, ageMax, lenMin, lenMax){

     samples      <- unique(dbeSimObj@ageStruc$estim$sample)
     times        <- as.character(unique(dbeSimObj@ageStruc$estim$time))
     spaces       <- as.character(unique(dbeSimObj@ageStruc$estim$space))
     technicals   <- as.character(unique(dbeSimObj@ageStruc$estim$technical))

    age <- dbeSimObj@ageStruc$estim
    len <- dbeSimObj@lenStruc$estim

    ages    <- ageMin:ageMax
    lens    <- lenMin:lenMax
    
    for(sm in samples){
        for(tm in times){
            for(sp in spaces){
                for(tc in technicals){
                     r1   <- subset(dbeSimObj@ageStruc$estim, sample == sm & time == tm & space == sp & technical == tc)
                     r2   <- subset(dbeSimObj@lenStruc$estim, sample == sm & time == tm & space == sp & technical == tc)
                     mis1 <- ages[which(!(ages %in% as.numeric(as.character(r1$age))))]
                     mis2 <- lens[which(!(lens %in% as.numeric(as.character(r2$len))))]
                     
                     age <- rbind(age, cbind(sample = rep(sm, length(mis1)), time = rep(tm, length(mis1)), space = rep(sp, length(mis1)),
                                             technical = rep(tc, length(mis1)), age = mis1, value = rep(0, length(mis1))))
                     len <- rbind(len, cbind(sample = rep(sm, length(mis2)), time = rep(tm, length(mis2)), space = rep(sp, length(mis2)),
                                             technical = rep(tc, length(mis2)), length = mis2, value = rep(0, length(mis2))))
                     
        }}}}
        
    age <- age[order(age$sample, age$time, age$space, age$technical, as.numeric(as.character(age$age))),]
    rownames(age) <- 1:dim(age)[1]
    age$value <- as.numeric(age$value)
    len <- len[order(len$sample, len$time, len$space, len$technical, as.numeric(as.character(len$length))),]
    rownames(len) <- 1:dim(len)[1]
    len$value <- as.numeric(len$value)
    
    dbeSimObj@ageStruc$estim <- age
    dbeSimObj@lenStruc$estim <- len
    
    return(dbeSimObj)}
    
    

fillGaps.dis <- function(dbeSimObj, ageMin, ageMax, lenMin, lenMax){

     samples      <- unique(dbeSimObj@lenStruc$estim$sample)
     times        <- as.character(unique(dbeSimObj@lenStruc$estim$time))
     spaces       <- as.character(unique(dbeSimObj@lenStruc$estim$space))
     technicals   <- as.character(unique(dbeSimObj@lenStruc$estim$technical))

    len <- dbeSimObj@lenStruc$estim

    lens    <- lenMin:lenMax
    
    for(sm in samples){
        for(tm in times){
            for(sp in spaces){
                for(tc in technicals){

                     r2   <- subset(dbeSimObj@lenStruc$estim, sample == sm & time == tm & space == sp & technical == tc)
         
                     mis2 <- lens[which(!(lens %in% as.numeric(as.character(r2$len))))]
                     

                     len <- rbind(len, cbind(sample = rep(sm, length(mis2)), time = rep(tm, length(mis2)), space = rep(sp, length(mis2)),
                                             technical = rep(tc, length(mis2)), length = mis2, value = rep(0, length(mis2))))
                     
        }}}}
        
    len <- len[order(len$sample, len$time, len$space, len$technical, as.numeric(as.character(len$length))),]
    rownames(len) <- 1:dim(len)[1]
    len$value <- as.numeric(len$value)
    

    dbeSimObj@lenStruc$estim <- len
    
    return(dbeSimObj)}
    
    