SortLoadings <-
function(Px, cutoff=min(apply(abs(Px),2,max))-1e-5){
  part <- Px
  ind <- 1
  
  for (r in 1:nrow(Px)){
    if (min(abs(Px[r,ind:ncol(Px)]))<cutoff){
      sorted <- sort(abs(part[r,]),decreasing=TRUE)
      IX <- colnames(sorted)
      Px[,ind:ncol(Px)] <- part[,IX]
      colnames(Px)[ind:ncol(Px)] <- IX
      indp <- min(which(abs(Px[r,ind:ncol(Px)])<cutoff))
      ind <- ind+indp-1
      part <- Px[,ind:ncol(Px)]
    }
    if ((r==nrow(Px))&(is.null(dim(part))==FALSE)){
      sorted <- sort(abs(part[r,]),decreasing=TRUE)
      IX <- colnames(sorted)
      Px[,ind:ncol(Px)] <- part[,IX]
      colnames(Px)[ind:ncol(Px)] <- IX
    }}
  return(Px)
}
