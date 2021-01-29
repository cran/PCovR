SortLoadings <-
  function(Px){
    Px_t = as.data.frame(t(Px))
    
    R <- ncol(Px_t)
    index <- 1
    
    maxi <- apply(abs(Px_t),1,which.max)
    sorted <- cbind(Px_t,maxi)
    order <- sorted[ order(abs(sorted[,R+1]),decreasing=F), ]
    IX <- rownames(order)
    sorted <- sorted[IX,]
    name <- NULL
    
    for (r in 1:ncol(Px_t)){
      part <- sorted[ order(abs(sorted[sorted[,R+1]==r, r]),decreasing=T) + index - 1, ]
      IX <- rownames(part)
      block <- sorted[IX,1:R]
      Px_t[index:(index+length(IX)-1),] <- block
      name <- c(name,IX)
      index <- index+length(IX)
    }
    rownames(Px_t) <- name
    Px = as.data.frame(t(Px_t))
    return(Px)
  }

