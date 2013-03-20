pcovr_est <- 
  function(X,Y,r,a,cross=FALSE,fold="LeaveOneOut"){
    library(MASS)
    library(ThreeWay)
    N <- nrow(X)
    
    #compute projector on X
    S <- t(X) %*% X
    if (det(S)<1e-12){ 
      l <- eigen(S)$values
      w <- which(l>1e-8*max(l))
      Sh <- eigen(S)$vectors[,w] %*% diag(l[w]^-1) %*% t(eigen(S)$vectors[,w])
      Hx <- X %*% Sh %*% t(X)
    } else {
      Hx <- X %*% solve(S) %*% t(X)
    }
    
    # Compute PCovR solution
    G <- a*X %*% t(X)*SUM(X)$ssq^-1+(1-a)*Hx %*% Y %*% t(Y) %*% Hx*SUM(Y)$ssq^-1
    G <- .5*(G+t(G)) #symmetrisation to avoid strange solutions
    Te <- eigen(G)$vectors[,1:r]
    W <- ginv(X) %*% Te
    Px <- t(Te) %*% X
    Py <- t(Te) %*% Y
    Rx2 <- SUM(Te %*% Px)$ssq*SUM(X)$ssq^-1
    Ry2 <- SUM(Te %*% Py)$ssq*SUM(Y)$ssq^-1
    B <- W %*% Py
    
    if (cross==TRUE){
      if (fold == "LeaveOneOut"){
        fold <- N
      }
      Yhatcv <- Y
      count <- 1
      LeaveOut <- rep(floor(N/fold), times = fold)
      LeaveOut[1] <- LeaveOut[1] + N%%fold
      CumulSum <- cumsum(LeaveOut)
      for (i in 1:fold){
        Xi <- X
        Xi <- Xi[-(count:CumulSum[i]),]
        Yi <- Y
        Yi <- array(Yi[-(count:CumulSum[i]),],c(nrow(X)-LeaveOut[i],ncol(Y)))
        
        S <- t(Xi) %*% Xi
        if (det(S)<1e-12){ 
          l <- eigen(S)$values
          w <- which(l>1e-8*max(l))
          Sh <- eigen(S)$vectors[,w] %*% diag(l[w]^-1) %*% t(eigen(S)$vectors[,w])
          Hx <- Xi %*% Sh %*% t(Xi)
        } else {
          Hx <- Xi %*% solve(S) %*% t(Xi)
        }
        
        G <- a*Xi %*% t(Xi)*SUM(Xi)$ssq^-1+(1-a)*Hx %*% Yi %*% t(Yi) %*% Hx*SUM(Yi)$ssq^-1
        G <- .5*(G+t(G))
        Ti <- eigen(G)$vectors[,1:r]
        Wi <- ginv(Xi) %*% Ti
        Bi <- Wi %*% t(Ti) %*% Yi
        Yhatcv[(count:CumulSum[i]),] <- X[(count:CumulSum[i]),] %*% Bi
        count <- count + LeaveOut[i]
      }
      Qy2 <- 1-SUM(Y-Yhatcv)$ssq*SUM(Y)$ssq^-1
      results <- list(W=W,B=B,Rx2=Rx2,Ry2=Ry2,Te=Te,Px=Px,Py=Py,Qy2=Qy2)
    } else {
      results <- list(W=W,B=B,Rx2=Rx2,Ry2=Ry2,Te=Te,Px=Px,Py=Py)
    }
    return(results)
  }