# function to find clusters
# written by Sarah Leavitt
# https://github.com/sarahleavitt/nbTransmission


clusterInfectors <- function(df, indIDVar, pVar,
                             clustMethod = c("n", "kd", "hc_absolute", "hc_relative"),
                             cutoff){
  
  df <- as.data.frame(df)
  #Creating variables with the individual ID
  indIDVar1 <- paste0(indIDVar, ".1")
  indIDVar2 <- paste0(indIDVar, ".2")
  
  #Checking that the named variables are in the data frame
  if(!indIDVar1 %in% names(df)){
    stop(paste0(indIDVar1, " is not in the data frame."))
  }
  if(!indIDVar2 %in% names(df)){
    stop(paste0(indIDVar2, " is not in the data frame."))
  }
  if(!pVar %in% names(df)){
    stop(paste0(pVar, " is not in the data frame."))
  }
  
  #Making sure clustMethod is correctly specified
  if(length(clustMethod) > 1){
    stop("Please provide a clustering method")
  }
  else if(!clustMethod %in% c("n", "kd", "hc_absolute", "hc_relative")){
    stop("clustMethod must be one of: n, kd, hc_absolute, hc_relative")
  }
  
  
  #Ranking the probabilities for each possible infector
  #Ties are set to the minimum rank of that group
  df <- df[order(df[, indIDVar2], -df[, pVar]), ]
  df$pRank <- stats::ave(-df[, pVar], df[, indIDVar2], 
                         FUN = function(x){
                           rank(x, ties.method = "min") 
                         })
  
  #Adding the number of infectors
  nInf <- as.data.frame(table(df[, indIDVar2]))
  names(nInf) <- c(indIDVar2, "nInfectors")
  df2 <- merge(df, nInf, by = indIDVar2, all = TRUE)
  
  #Splitting the data frame to those who have one infector and multiple infectors
  multInf <- df2[df2$nInfectors > 1, ]
  oneInf <- df2[df2$nInfectors == 1, ]
  
  ## Using a constant number of infectors ##
  if(clustMethod == "n"){
    
    clustRes <- multInf
    clustRes$cluster <- ifelse(clustRes$pRank <= cutoff, 1, 2)
    
  }
  
  ## Using kernel density estimation ##
  if(clustMethod == "kd"){
    
    clustRes1 <- dplyr::group_by(multInf, !!rlang::sym(indIDVar2))
    clustRes <- dplyr::group_modify(clustRes1, ~ findClustersKD(.x, pVar = pVar,
                                                                cutoff = cutoff))
    
  }
  
  ## Using hierarchical clustering ##
  if(grepl("hc", clustMethod)){
    
    clustRes1 <- dplyr::group_by(multInf, !!rlang::sym(indIDVar2))
    clustRes <- dplyr::group_modify(clustRes1, ~ findClustersHC(.x, pVar,
                                                                cutoff = cutoff,
                                                                clustMethod = clustMethod))
    
  }
  
  #Removing tibble formatting
  clustRes <- as.data.frame(dplyr::ungroup(clustRes))
  
  #Combining the clustering for those with more than one infector with those
  #who have one infector
  if(nrow(oneInf) > 0){
    #If there is one infector, it should be in the top cluster
    oneInf$cluster <- 1
    clustRes2 <- rbind(oneInf, clustRes)
  }else{
    clustRes2 <- clustRes
  }
  
  #Making the clusters a factor variable
  clustRes2$cluster <- factor(clustRes2$cluster, levels = c(1, 2))
  #Removing nInfectors variable
  clustRes2$nInfectors <- NULL
  
  return(clustRes2)
}



