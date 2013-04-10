SortLoadings <- 
  function(Px, cutoff = .20){
  part <- Px
  ind <- 1
  
  for (r in 1:nrow(Px)){
    sorted <- sort(abs(part[r,]),decreasing=TRUE)
    IX <- colnames(sorted)
    Px[,ind:ncol(Px)] <- part[,IX]
    colnames(Px)[ind:ncol(Px)] <- IX
    indp <- min(which(abs(Px[r,ind:ncol(Px)])<cutoff))
    ind <- ind+indp-1
    part <- Px[,ind:ncol(Px)]
  }
  return(Px)
}